# File: backend/app/core/assembly/hierarchical_assembler.py
# Version: v0.1.4

"""
Hierarchical assembler: recursively subdivide a target sequence into fragments,
compute overlaps between each layer’s children via GA, and build a tree of
FragmentNode that you can later traverse or export via separate utilities.
"""

from dataclasses import dataclass, field
from pathlib import Path
import math
from typing import List, Tuple, Dict, Optional

from backend.app.core.config import LevelConfig, load_levels_config
from backend.app.core.assembly.ga_overlap_selector import GAOverlapSelector


@dataclass
class FragmentNode:
    """
    Represents one node in the assembly tree.
    """
    level:       int
    start:       int                # absolute start in the root sequence
    end:         int                # absolute end in the root sequence
    seq:         str                # slice of the full root sequence [start:end]
    overlaps:    List[Tuple[str,int,int]] = field(default_factory=list)
    ga_log:      List[float]            = field(default_factory=list)
    children:    List["FragmentNode"]   = field(default_factory=list)


class HierarchicalAssembler:
    """
    Drives the recursive fragment-split → GA-overlap → recurse process.
    """
    def __init__(self, levels_cfg: Dict[int, LevelConfig]):
        """
        Args:
            levels_cfg: mapping level index → LevelConfig (loaded from JSON).
        """
        self.levels = levels_cfg

    def build(self, full_seq: str) -> FragmentNode:
        """
        Start the recursion at level 0, covering the entire sequence.
        """
        return self._build_node(full_seq, level=0, start=0, end=len(full_seq))

    def _build_node(
        self,
        full_seq: str,
        level: int,
        start: int,
        end: int
    ) -> FragmentNode:
        """
        Recursively split `full_seq[start:end]` at this `level`,
        run GA to compute overlaps, then recurse on children using
        the full fragment (body + overlaps) for each.
        """
        seq = full_seq[start:end]
        cfg = self.levels.get(level)
        node = FragmentNode(level=level, start=start, end=end, seq=seq)

        # Stop if no config for this level or fragment too small
        if cfg is None or (end - start) < cfg.fragment_min_size:
            return node

        # 1) Determine how many children to split into
        L       = end - start
        ideal_n = max(1, round(L / cfg.fragment_max_size))
        n       = min(max(ideal_n, cfg.min_children), cfg.max_children)

        # 2) Compute equally spaced child‐body boundaries
        bounds: List[Tuple[int,int]] = []
        for i in range(n):
            s_i = start + math.floor(i   * L / n)
            e_i = start + math.floor((i+1) * L / n)
            bounds.append((s_i, e_i))

        # 3) Run GAOverlapSelector on the raw bodies
        bodies    = [full_seq[s:e] for s,e in bounds]
        positions = bounds
        ga = GAOverlapSelector(
            full_sequence    = full_seq,
            oligo_positions  = positions,
            oligo_seqs       = bodies,
            overlap_min_size = cfg.overlap_min_size,
            overlap_max_size = cfg.overlap_max_size
        )
        print(f"Level {level}: splitting into {n} children, running GA…")
        best, log = ga.evolve()
        node.overlaps = best.overlaps
        node.ga_log   = log

        # 4) Compute each child's **extended** bounds (body + its overlaps)
        extended_bounds: List[Tuple[int,int]] = []
        for idx, (s_i, e_i) in enumerate(bounds):
            # upstream overlap for child idx
            prev_ov = best.overlaps[idx-1] if idx > 0 else None
            # downstream overlap for child idx
            next_ov = best.overlaps[idx]   if idx < len(bounds)-1 else None

            new_start = prev_ov[1] if prev_ov else s_i
            new_end   = next_ov[2] if next_ov else e_i
            extended_bounds.append((new_start, new_end))

        # 5) Recurse on each **extended** fragment
        for s_i, e_i in extended_bounds:
            child = self._build_node(full_seq, level+1, s_i, e_i)
            node.children.append(child)

        return node



# -----------------------------------------------------------------------------
# Example usage:
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # 1) load hierarchical levels from JSON:
    levels = load_levels_config(Path("app/config/levels.json"))

    # 2) read one sequence from FASTA:
    from Bio import SeqIO
    rec    = next(SeqIO.parse("data/input_sequences_rnd_9000.fasta","fasta"))
    seq    = str(rec.seq).upper()
    root_id = rec.id

    # 3) build the hierarchy:
    asm  = HierarchicalAssembler(levels)
    tree = asm.build(seq)

    # 4) export everything via export_utils:
    from backend.app.core.assembly.export_utils import (
        export_fragments_to_fasta,
        export_tree_to_json,
        export_tree_to_html,
        export_clusters_and_fragments_html
    )

    export_fragments_to_fasta(tree, Path("data/tree_fragments.fasta"), root_id)
    export_tree_to_json(tree,      Path("data/tree_fragments.json"),  root_id)
    export_tree_to_html(tree,      Path("data/tree_fragments.html"),  root_id)
    export_clusters_and_fragments_html(
        tree,
        Path("data/tree_clusters_fragments.html"),
        root_id
    )
