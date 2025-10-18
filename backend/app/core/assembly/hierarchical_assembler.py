# File: backend/app/core/assembly/hierarchical_assembler.py
# Version: v0.8.1
"""
Hierarchical assembler: recursively subdivide the target into fragments,
optimize overlaps between children via GA, and build a FragmentNode tree.

v0.8.1 — User-facing base-1 labels
- Child fragment_id now uses P-<parent+1>, F-<idx+1> (was base-0).
- Log messages mentioning "F-<n>" are printed base-1.

v0.8.0 — Enforce per-junction minimum overlaps from the planner + auto-fix loop
- The end-aware planner now returns (positions, hat_min_overlaps).
- Pass those minimums into GA (min_overlap_per_junction) so GA cannot shrink
  children below the planner’s expectations.
- After GA, validate FULL child lengths. If any child violates [min..max],
  increase the adjacent junction-minimums just enough to fix the deficit and
  re-run GA (bounded retries). This prevents <min children like 91 bp.
"""

from __future__ import annotations

import logging
from typing import List, Tuple, Dict, Optional

from backend.app.core.models.fragment_node import FragmentNode
from backend.app.config.config_level import LevelConfig
from backend.app.config.config_global import GlobalConfig
from backend.app.core.assembly.ga_overlap_selector import (
    GAOverlapSelector,
    NoOverlapCandidatesError,
)
from backend.app.core.assembly.endaware_planner import plan_bodies_endaware_with_prescan


class AssemblerError(RuntimeError):
    pass


MAX_GA_RETRIES = 12
MAX_FIX_ROUNDS = 6  # times we'll bump per-junction mins and retry GA


class HierarchicalAssembler:
    def __init__(
        self,
        levels_cfg: Dict[int, LevelConfig],
        global_cfg: Optional[GlobalConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        self.levels = levels_cfg
        self.global_cfg = global_cfg or GlobalConfig()
        self.sequence_id: Optional[str] = None
        self.log = logger or logging.getLogger(__name__)

    # --------------------- Public API ---------------------

    def build(self, full_seq: str, root_id: str) -> FragmentNode:
        """Start recursion at level 0, covering the entire sequence."""
        self.sequence_id = root_id

        l0 = self.levels.get(0)
        if l0 is None:
            raise ValueError("levels.json must define level=0")

        L0 = len(full_seq)
        if not (l0.fragment_min_size <= L0 <= l0.fragment_max_size):
            raise ValueError(
                f"[L0] input sequence length {L0} out of fragment range "
                f"({l0.fragment_min_size}..{l0.fragment_max_size})"
            )

        root = FragmentNode(
            fragment_id=f"S-{root_id}_L-0_P-0_F-0",
            level=0,
            start=0,
            end=len(full_seq),
            seq=full_seq,
            strand="+",
            overlap_prev="",
            overlap_next="",
            is_oligo=False,
            ga_log=[],
            children=[],
        )

        root.children = self._build_children(
            full_seq=full_seq, level=0, start=0, end=len(full_seq),
            parent_frag_index=0, parent_node=root
        )
        return root

    # --------------------- Recursion ---------------------

    def _build_children(
        self,
        full_seq: str,
        level: int,
        start: int,
        end: int,
        parent_frag_index: int,
        parent_node: FragmentNode,
    ) -> List[FragmentNode]:
        cfg = self.levels.get(level)
        L = end - start

        if cfg is None:
            return []

        if not (cfg.fragment_min_size <= L <= cfg.fragment_max_size):
            raise ValueError(
                f"[L{level}] fragment length {L} out of range "
                f"({cfg.fragment_min_size}..{cfg.fragment_max_size})"
            )

        # --- End-aware, pre-scan (positions + per-junction **min** overlaps) ---
        positions, min_ov_lens = plan_bodies_endaware_with_prescan(
            full_seq=full_seq,
            start=start,
            end=end,
            cfg=cfg,
            tm_method=self.global_cfg.tm_method,
            tm_params=self.global_cfg.tm_params_dict(),
        )

        def make_ga(curr_positions: List[Tuple[int, int]],
                    min_j: Optional[List[int]]) -> GAOverlapSelector:
            ga = GAOverlapSelector(
                full_sequence=full_seq,
                oligo_positions=curr_positions,
                oligo_seqs=[full_seq[s:e] for s, e in curr_positions],
                overlap_min_size=cfg.overlap_min_size,
                overlap_max_size=cfg.overlap_max_size,
                tm_min=cfg.overlap_tm_min,
                tm_max=cfg.overlap_tm_max,
                gc_min=cfg.overlap_gc_min,
                gc_max=cfg.overlap_gc_max,
                run_min=cfg.overlap_run_min,
                run_max=cfg.overlap_run_max,
                disallowed_motifs=cfg.overlap_disallowed_motifs,
                allowed_motifs=cfg.overlap_allowed_motifs,
                max_gap_allowed=getattr(cfg, "max_gap_allowed", 0),
                enforce_span=True,  # must straddle the junction
                tm_method=self.global_cfg.tm_method,
                tm_params=self.global_cfg.tm_params_dict(),
                ga_params=self.global_cfg.ga_params(),
                min_child_full_len=cfg.min_children_size,
                max_child_full_len=cfg.max_children_size,
                # CRITICAL: prevent GA from picking overlaps shorter than planner’s minima
                min_overlap_per_junction=min_j,
            )
            ga.set_cpu_fraction(self.global_cfg.cpu_workers_fraction)
            return ga

        attempts = 0
        last_err: Optional[NoOverlapCandidatesError] = None

        # Keep a working copy of minima we can bump if needed
        work_min = list(min_ov_lens)

        while attempts < MAX_GA_RETRIES:
            try:
                self.log.info(
                    "Level %d: GA on %d junction(s); overlap %d-%d",
                    level, len(positions) - 1, cfg.overlap_min_size, cfg.overlap_max_size,
                )
                # Inner “fix-up” loop: if after GA any child violates size bounds, bump mins and retry GA
                fix_round = 0
                while fix_round <= MAX_FIX_ROUNDS:
                    ga = make_ga(positions, work_min)
                    best, progress_log = ga.evolve()

                    # Attach GA progress to the **parent/cluster** node (for reporting)
                    parent_node.ga_log = list(progress_log)
                    parent_node.ga_detail = list(ga.progress_detail)
                    parent_node.ga_cluster_id = f"S-{self.sequence_id}_L-{level}_S-{start}_E-{end}"

                    # Validate FULL child lengths against cfg bounds
                    bad_idx = self._find_first_child_size_violation(full_seq, positions, best.overlaps, cfg)
                    if bad_idx is None:
                        # All good, build children nodes
                        children = self._materialize_children(full_seq, positions, best.overlaps, level, parent_frag_index)
                        return children
                    else:
                        # Too short/long → bump adjacent per-junction minimum overlaps and retry GA
                        i, length, vtype = bad_idx
                        deficit = (cfg.min_children_size - length) if vtype == "too_short" else (length - cfg.max_children_size)
                        self._bump_minima_for_child(i, deficit, work_min, cfg)
                        # NOTE: print F-<idx+1> (base-1)
                        self.log.warning(
                            "[L%d] Child F-%d length %d is %s bound [%d..%d] ⇒ bump mins → %s; fix_round %d/%d",
                            level, i + 1, length, "below" if vtype == "too_short" else "above",
                            cfg.min_children_size, cfg.max_children_size, work_min, fix_round+1, MAX_FIX_ROUNDS
                        )
                        fix_round += 1

                # If we drop out of fix loop, give up with a helpful error
                raise AssemblerError(
                    f"[L{level}] Unable to satisfy child size bounds after {MAX_FIX_ROUNDS} fix rounds; "
                    f"consider relaxing overlap or child size constraints."
                )

            except NoOverlapCandidatesError as err:
                attempts += 1
                last_err = err
                self.log.warning(
                    ("Level %d: zero candidates at junction %d (cut %d|%d). "
                     "Reasons: %s — retry %d/%d"),
                    level, err.junction_index, err.left_end, err.right_start,
                    err.formatted_reasons(), attempts, MAX_GA_RETRIES,
                )
                positions = self._jitter_positions(positions, radius=min(attempts, 5))
                # Also relax local minima slightly near the reported junction, if we have them
                if work_min and 0 <= err.junction_index < len(work_min):
                    work_min[err.junction_index] = max(cfg.overlap_min_size, work_min[err.junction_index] - 1)
        else:
            raise AssemblerError(self._format_infeasible_message(level, last_err, cfg))

    # --------------------- Helpers ---------------------

    def _materialize_children(
        self,
        full_seq: str,
        positions: List[Tuple[int, int]],
        overlaps: List[Tuple[str, int, int]],
        level: int,
        parent_frag_index: int,
    ) -> List[FragmentNode]:
        n = len(positions)
        children: List[FragmentNode] = []
        for idx, (s_i, e_i) in enumerate(positions):
            prev_ov = overlaps[idx - 1] if idx > 0 else None
            next_ov = overlaps[idx] if idx < n - 1 else None

            ext_s = prev_ov[1] if prev_ov else s_i
            ext_e = next_ov[2] if next_ov else e_i
            frag_seq = full_seq[ext_s:ext_e]
            strand = "+" if (idx % 2 == 0) else "-"

            prev_seq = full_seq[prev_ov[1]:prev_ov[2]] if prev_ov else ""
            next_seq = full_seq[next_ov[1]:next_ov[2]] if next_ov else ""

            is_leaf = (self.levels.get(level + 1) is None)

            # USER-FACING base-1 labels for P and F:
            fragment_id = f"S-{self.sequence_id}_L-{level+1}_P-{parent_frag_index+1}_F-{idx+1}"

            node = FragmentNode(
                fragment_id=fragment_id,
                level=level + 1,
                start=ext_s,
                end=ext_e,
                seq=frag_seq,
                strand=strand,
                overlap_prev=prev_seq,
                overlap_next=next_seq,
                is_oligo=is_leaf,
                ga_log=[],
                children=[],
            )

            # Recurse: child becomes parent for next level (parent_frag_index remains numeric)
            node.children = self._build_children(
                full_seq, level + 1, ext_s, ext_e, parent_frag_index=idx, parent_node=node
            )
            children.append(node)
        return children

    def _find_first_child_size_violation(
        self,
        full_seq: str,
        positions: List[Tuple[int, int]],
        overlaps: List[Tuple[str, int, int]],
        cfg: LevelConfig,
    ) -> Optional[Tuple[int, int, str]]:
        """
        Return (index, length, 'too_short'|'too_long') of the first violating child,
        or None if all are within [min..max].
        """
        n = len(positions)
        for idx, (s_i, e_i) in enumerate(positions):
            prev_ov = overlaps[idx - 1] if idx > 0 else None
            next_ov = overlaps[idx] if idx < n - 1 else None
            ext_s = prev_ov[1] if prev_ov else s_i
            ext_e = next_ov[2] if next_ov else e_i
            Lfull = ext_e - ext_s
            if Lfull < cfg.min_children_size:
                return (idx, Lfull, "too_short")
            if Lfull > cfg.max_children_size:
                return (idx, Lfull, "too_long")
        return None

    def _bump_minima_for_child(
        self,
        child_idx: int,
        deficit: int,
        work_min: List[int],
        cfg: LevelConfig,
    ) -> None:
        """
        Increase the minimum overlap(s) adjacent to child `child_idx` just enough to
        fix the size violation by spreading the delta over the two sides.
        """
        if deficit <= 0:
            return
        nJ = len(work_min)
        # Distribute roughly evenly to left/right junctions
        left_add = (deficit + 1) // 2
        right_add = deficit // 2
        if child_idx > 0:
            jL = child_idx - 1
            work_min[jL] = min(cfg.overlap_max_size, max(cfg.overlap_min_size, work_min[jL] + left_add))
        if child_idx < nJ:
            jR = child_idx
            work_min[jR] = min(cfg.overlap_max_size, max(cfg.overlap_min_size, work_min[jR] + right_add))

    def _jitter_positions(
        self,
        positions: List[Tuple[int, int]],
        radius: int = 1,
    ) -> List[Tuple[int, int]]:
        if radius <= 0 or len(positions) <= 2:
            return positions
        new_pos = positions[:]
        for j in range(len(positions) - 1):
            ls, le = new_pos[j]
            rs, re = new_pos[j + 1]
            cut = le
            moved = False
            for d in range(1, radius + 1):
                for delta in (d, -d):
                    nc = cut + delta
                    if ls < nc < re:
                        new_pos[j] = (ls, nc)
                        new_pos[j + 1] = (nc, re)
                        moved = True
                        break
                if moved:
                    break
        return new_pos

    def _format_infeasible_message(
        self,
        level: int,
        err: Optional[NoOverlapCandidatesError],
        cfg: LevelConfig,
    ) -> str:
        if err is None:
            return (f"[L{level}] No feasible overlap configuration after {MAX_GA_RETRIES} retries; "
                    "consider relaxing overlap constraints.")
        reasons = getattr(err, "reasons", {}) or {}
        reason_str = ", ".join(f"{k}:{v}" for k, v in reasons.items()) if reasons else "no-counters"
        sugg = self._suggestions_from_reasons(reasons, cfg)
        hints = " ".join(sugg) if sugg else "Consider relaxing Tm/GC/run/motif or widening overlap length range."
        return (
            f"[L{level}] Infeasible overlaps at junction {err.junction_index} (cut {err.left_end}|{err.right_start}) "
            f"after {MAX_GA_RETRIES} retries. Reasons: {reason_str}. {hints}"
        )

    @staticmethod
    def _suggestions_from_reasons(reasons: Dict[str, int], cfg: LevelConfig) -> List[str]:
        if not reasons:
            return []
        tips: List[str] = []
        total = sum(reasons.values()) or 1

        def dominant(key: str, frac: float = 0.25) -> bool:
            return reasons.get(key, 0) / total >= frac

        if dominant("tm_low"):
            tips.append(f"Tip: lower overlap_tm_min (currently {cfg.overlap_tm_min}) by ~3–5 °C.")
        if dominant("tm_high"):
            tips.append(f"Tip: raise overlap_tm_max (currently {cfg.overlap_tm_max}) by ~2–3 °C.")
        if dominant("gc_low"):
            tips.append(f"Tip: lower overlap_gc_min (currently {cfg.overlap_gc_min}) by ~3–5%%.")
        if dominant("gc_high"):
            tips.append(f"Tip: raise overlap_gc_max (currently {cfg.overlap_gc_max}) by ~3–5%%.")
        if dominant("run_high"):
            tips.append(f"Tip: relax overlap_run_max (currently {cfg.overlap_run_max}).")
        if dominant("run_low"):
            tips.append("Tip: if you set overlap_run_min, consider lowering it.")
        if dominant("motif_disallowed"):
            tips.append("Tip: review overlap_disallowed_motifs; remove overly generic motifs.")
        if dominant("length_window"):
            tips.append(
                f"Tip: widen overlap length range (currently {cfg.overlap_min_size}–{cfg.overlap_max_size} nt) "
                "or allow a higher body count at this level."
            )
        if not tips:
            tips.append("Tip: relax one of Tm/GC/run/motif constraints, or widen overlap length range.")
        return tips
