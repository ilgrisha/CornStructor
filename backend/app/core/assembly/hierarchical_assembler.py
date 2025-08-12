# File: backend/app/core/assembly/hierarchical_assembler.py
# Version: v0.4.2

"""
Hierarchical assembler: recursively subdivide a target sequence into fragments,
compute overlaps between each layer’s children via GA, and build a tree of
FragmentNode that you can later traverse or export via separate utilities.

v0.4.2:
- Fix regression where _build_children (and helpers) were omitted in v0.4.1.
- Keeps import switch to use export_all_levels (per-level report set).
- Integrates division-point retries, per-level gap, and advanced overlap filters.
"""

from pathlib import Path
import math
import logging
from typing import List, Tuple, Dict, Optional

from backend.app.core.models.fragment_node import FragmentNode
from backend.app.core.config import LevelConfig, load_levels_config
from backend.app.core.config_global import GlobalConfig, load_global_config
from backend.app.core.assembly.ga_overlap_selector import (
    GAOverlapSelector,
    NoOverlapCandidatesError,
)
from backend.app.core.visualization.oligo_html_visualizer import reverse_complement

from backend.app.core.export.fasta_exporter import export_fragments_to_fasta
from backend.app.core.export.json_exporter import export_tree_to_json
from backend.app.core.visualization.tree_html_exporter import export_tree_to_html
from backend.app.core.visualization.cluster_html_report import export_all_levels


class AssemblerError(RuntimeError):
    """Raised when hierarchical assembly cannot satisfy overlap constraints after retries."""


class HierarchicalAssembler:
    """
    Drives the recursive fragment-split → GA-overlap → recurse process.
    """

    def __init__(
        self,
        levels_cfg: Dict[int, LevelConfig],
        global_cfg: Optional[GlobalConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """
        Args:
            levels_cfg: mapping level index → LevelConfig (loaded from JSON).
            global_cfg: global parameters (GA knobs, division-point adjustment, Tm params).
            logger:     optional logger; if None, a module logger is used.
        """
        self.levels = levels_cfg
        self.global_cfg = global_cfg or GlobalConfig()
        self.sequence_id: Optional[str] = None
        self.log = logger or logging.getLogger(__name__)

    # --------------------- Public API ---------------------

    def build(self, full_seq: str, root_id: str) -> FragmentNode:
        """
        Start the recursion at level 0, covering the entire sequence.
        The root node is always sense strand ('+').
        """
        self.sequence_id = root_id

        root_node = FragmentNode(
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

        root_node.children = self._build_children(
            full_seq, level=0, start=0, end=len(full_seq), parent_frag_index=0
        )

        return root_node

    # --------------------- Internals ---------------------

    def _run_ga_or_adjust(
        self,
        full_seq: str,
        level: int,
        positions: List[Tuple[int, int]],
        bodies: List[str],
        cfg: LevelConfig,
    ):
        """
        Try GA; if filters eliminate all candidates at a junction, attempt to
        move that division cut several times per global config until success
        or exhaustion, while respecting next-level child length bounds.

        Returns:
            (best_chromosome, progress_log, adjusted_positions, adjusted_bodies)
        """
        next_cfg = self.levels.get(level + 1)  # bounds for children (if any)
        attempts_left = self.global_cfg.max_total_adjustment_rounds

        def make_ga(curr_positions: List[Tuple[int, int]]) -> GAOverlapSelector:
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
                max_gap_allowed=cfg.max_gap_allowed,
                enforce_span=self.global_cfg.enforce_span_across_junction,
                tm_method=self.global_cfg.tm_method,
                tm_params=self.global_cfg.tm_params_dict(),
                ga_params=self.global_cfg.ga_params(),
            )
            ga.set_cpu_fraction(self.global_cfg.cpu_workers_fraction)
            return ga

        while True:
            try:
                ga = make_ga(positions)
                self.log.info(
                    "Level %d: GA on %d junction(s); overlap %d-%d, Tm %s, GC %s, run %s, gap<=%d",
                    level,
                    len(positions) - 1,
                    cfg.overlap_min_size,
                    cfg.overlap_max_size,
                    f"{cfg.overlap_tm_min}-{cfg.overlap_tm_max}",
                    f"{cfg.overlap_gc_min}-{cfg.overlap_gc_max}",
                    f"{cfg.overlap_run_min}-{cfg.overlap_run_max}",
                    cfg.max_gap_allowed or 0,
                )
                best, log = ga.evolve()
                return best, log, positions, [full_seq[s:e] for s, e in positions]

            except NoOverlapCandidatesError as err:
                # Log failing junction and constraint counts
                self.log.warning(
                    ("Level %d: zero candidates at junction %d (cut %d|%d). "
                     "Reasons: %s"),
                    level,
                    err.junction_index,
                    err.left_end,
                    err.right_start,
                    err.formatted_reasons(),
                )
                if attempts_left <= 0:
                    raise AssemblerError(
                        f"Level {level}: cannot satisfy overlap constraints at cut "
                        f"{err.left_end}|{err.right_start} after "
                        f"{self.global_cfg.max_total_adjustment_rounds} adjustment round(s). "
                        f"Reasons: {err.formatted_reasons()}"
                    ) from err

                # Try moving the specific division point while keeping children within size bounds
                j = err.junction_index
                positions2 = positions[:]
                cut = positions2[j][1]  # == positions2[j+1][0]
                deltas = []
                step = self.global_cfg.division_adjust_step_bp
                for k in range(1, self.global_cfg.division_adjust_attempts_per_junction + 1):
                    mag = step * k
                    deltas.extend([+mag, -mag])

                for delta in deltas:
                    new_cut = cut + delta
                    left_s, left_e = positions2[j]
                    right_s, right_e = positions2[j + 1]
                    if not (left_s < new_cut < right_e):
                        continue

                    if next_cfg:
                        left_len = new_cut - left_s
                        right_len = right_e - new_cut
                        if left_len < next_cfg.fragment_min_size or left_len > next_cfg.fragment_max_size:
                            continue
                        if right_len < next_cfg.fragment_min_size or right_len > next_cfg.fragment_max_size:
                            continue

                    positions_try = positions2[:]
                    positions_try[j] = (left_s, new_cut)
                    positions_try[j + 1] = (new_cut, right_e)

                    try:
                        ga_try = make_ga(positions_try)
                        self.log.info(
                            "Level %d: retry GA after shifting junction %d by %+d bp → cut @ %d",
                            level,
                            j,
                            delta,
                            new_cut,
                        )
                        best, log = ga_try.evolve()
                        return best, log, positions_try, [full_seq[s:e] for s, e in positions_try]
                    except NoOverlapCandidatesError as err2:
                        self.log.debug(
                            "Level %d: still failing after delta %+d at junction %d; reasons: %s",
                            level,
                            delta,
                            j,
                            err2.formatted_reasons(),
                        )
                        continue

                attempts_left -= 1
                positions = self._jitter_all_cuts(
                    positions, next_cfg, self.global_cfg.division_adjust_step_bp // 2
                )

    def _jitter_all_cuts(
        self,
        positions: List[Tuple[int, int]],
        next_cfg: Optional[LevelConfig],
        max_shift: int,
    ) -> List[Tuple[int, int]]:
        """
        Apply a small left/right jitter to each internal cut to escape synchronized dead-ends.
        Always respects child length bounds if next_cfg is present.
        """
        if max_shift <= 0 or len(positions) <= 2:
            return positions

        new_pos = positions[:]
        for j in range(len(positions) - 1):
            left_s, left_e = new_pos[j]
            right_s, right_e = new_pos[j + 1]
            cut = left_e
            for delta in (max_shift, -max_shift):
                new_cut = cut + delta
                if not (left_s < new_cut < right_e):
                    continue
                if next_cfg:
                    left_len = new_cut - left_s
                    right_len = right_e - new_cut
                    if not (next_cfg.fragment_min_size <= left_len <= next_cfg.fragment_max_size):
                        continue
                    if not (next_cfg.fragment_min_size <= right_len <= next_cfg.fragment_max_size):
                        continue
                new_pos[j] = (left_s, new_cut)
                new_pos[j + 1] = (new_cut, right_e)
                break
        return new_pos

    def _build_children(
        self,
        full_seq: str,
        level: int,
        start: int,
        end: int,
        parent_frag_index: int,
    ) -> List[FragmentNode]:
        """
        Helper to build child nodes of fragment [start:end] at given level.
        Returns a list of FragmentNode for this level’s children.
        """
        cfg = self.levels.get(level)
        L = end - start

        # If no further subdivision, this is a leaf oligo
        if cfg is None or L < cfg.fragment_min_size:
            return []

        # 1) determine number of children
        ideal_n = max(1, round(L / cfg.fragment_max_size))
        n = min(max(ideal_n, cfg.min_children), cfg.max_children)

        # 2) Compute equally spaced child‐body boundaries
        bounds: List[Tuple[int, int]] = []
        for i in range(n):
            s_i = start + math.floor(i * L / n)
            e_i = start + math.floor((i + 1) * L / n)
            bounds.append((s_i, e_i))

        bodies = [full_seq[s:e] for s, e in bounds]
        positions = bounds

        # 3) Run GA (with possible division-point adjustments on failure)
        best, log, positions, bodies = self._run_ga_or_adjust(
            full_seq=full_seq, level=level, positions=positions, bodies=bodies, cfg=cfg
        )

        # 4) compute extended bounds + overlap metadata per child
        children: List[FragmentNode] = []
        n = len(positions)
        for idx, ((s_i, e_i), body_seq) in enumerate(zip(positions, bodies)):
            prev_ov = best.overlaps[idx - 1] if idx > 0 else None
            next_ov = best.overlaps[idx] if idx < n - 1 else None

            ext_s = prev_ov[1] if prev_ov else s_i
            ext_e = next_ov[2] if next_ov else e_i
            frag_seq = full_seq[ext_s:ext_e]

            strand = "+" if (idx % 2 == 0) else "-"

            prev_seq = prev_ov[0] if prev_ov else ""
            next_seq = next_ov[0] if next_ov else ""
            if strand == "-":
                prev_seq = reverse_complement(prev_seq)
                next_seq = reverse_complement(next_seq)

            next_cfg = self.levels.get(level + 1)
            is_leaf = next_cfg is None or (ext_e - ext_s) < next_cfg.fragment_min_size

            fragment_id = f"S-{self.sequence_id}_L-{level+1}_P-{parent_frag_index}_F-{idx}"

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
                ga_log=log,
                children=[],
            )

            node.children = self._build_children(
                full_seq, level + 1, ext_s, ext_e, parent_frag_index=idx
            )
            children.append(node)

        return children


# -----------------------------------------------------------------------------
# Example usage (manual test)
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # 1) load hierarchical levels from JSON:
    levels = load_levels_config(Path("app/config/levels.json"))
    gconf = load_global_config(Path("app/config/globals.json"))

    # 2) read one sequence from FASTA:
    from Bio import SeqIO

    rec = next(SeqIO.parse("backend/data/input_sequences_rnd_9000.fasta", "fasta"))
    seq = str(rec.seq).upper()
    root_id = rec.id

    # 3) build the hierarchy:
    asm = HierarchicalAssembler(levels, gconf)
    root = asm.build(seq, root_id)

    # 4) export everything:
    export_fragments_to_fasta(root, Path("backend/data/out/tree_fragments.fasta"), root_id)
    export_tree_to_json(root, Path("backend/data/out/tree_fragments.json"), root_id)
    export_tree_to_html(root, Path("backend/data/out/tree_fragments.html"), root_id)
    export_all_levels(root, Path("backend/data/out/cluster_reports"), levels, gconf)
