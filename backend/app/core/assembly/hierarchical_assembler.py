# File: backend/app/core/assembly/hierarchical_assembler.py
# Version: v0.2.0

"""
Hierarchical assembler: recursively subdivide a target sequence into fragments,
compute overlaps between each layer’s children via GA, and build a tree of
FragmentNode that you can later traverse or export via separate utilities.
"""

from dataclasses import dataclass, field
from pathlib import Path
import math
from typing import List, Tuple, Dict, Optional

from backend.app.core.models.fragment_node import FragmentNode
from backend.app.core.config import LevelConfig, load_levels_config
from backend.app.core.assembly.ga_overlap_selector import GAOverlapSelector
from backend.app.core.visualization.oligo_html_visualizer import reverse_complement

from backend.app.core.export.fasta_exporter import export_fragments_to_fasta
from backend.app.core.export.json_exporter import export_tree_to_json
from backend.app.core.visualization.tree_html_exporter import export_tree_to_html
from backend.app.core.visualization.cluster_html_report import export_clusters_and_fragments_html


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
        self.sequence_id = None  # will be set during build

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
            strand='+',
            overlap_prev="",
            overlap_next="",
            is_oligo=False,
            ga_log=[],
            children=[]
        )

        root_node.children = self._build_children(
            full_seq,
            level=0,
            start=0,
            end=len(full_seq),
            parent_frag_index=0
        )

        return root_node

    def _build_children(
        self,
        full_seq: str,
        level: int,
        start: int,
        end: int,
        parent_frag_index: int
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

        # 3) run GA to pick overlaps
        ga = GAOverlapSelector(
            full_sequence=full_seq,
            oligo_positions=positions,
            oligo_seqs=bodies,
            overlap_min_size=cfg.overlap_min_size,
            overlap_max_size=cfg.overlap_max_size
        )
        print(f"Level {level}: splitting into {n} fragments, running GA…")
        best, log = ga.evolve()

        # 4) compute extended bounds + overlap metadata per child
        children: List[FragmentNode] = []
        for idx, (s_i, e_i) in enumerate(bounds):
            # upstream/downstream overlaps (raw)
            prev_ov = best.overlaps[idx - 1] if idx > 0 else None
            next_ov = best.overlaps[idx] if idx < n - 1 else None

            # calculate full fragment bounds including overlaps
            ext_s = prev_ov[1] if prev_ov else s_i
            ext_e = next_ov[2] if next_ov else e_i
            frag_seq = full_seq[ext_s:ext_e]

            # assign strand alternation per child index
            strand = '+' if (idx % 2 == 0) else '-'

            # orient overlap sequences to this strand
            prev_seq = prev_ov[0] if prev_ov else ""
            next_seq = next_ov[0] if next_ov else ""
            if strand == '-':
                prev_seq = reverse_complement(prev_seq)
                next_seq = reverse_complement(next_seq)

            # determine if this fragment should be further subdivided
            next_cfg = self.levels.get(level + 1)
            is_leaf = next_cfg is None or (ext_e - ext_s) < next_cfg.fragment_min_size

            # construct fragment_id
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
                children=[]
            )

            # recurse to build its children
            node.children = self._build_children(
                full_seq, level + 1, ext_s, ext_e, parent_frag_index=idx
            )
            children.append(node)

        return children


# -----------------------------------------------------------------------------
# Example usage:
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # 1) load hierarchical levels from JSON:
    levels = load_levels_config(Path("app/config/levels.json"))

    # 2) read one sequence from FASTA:
    from Bio import SeqIO
    rec = next(SeqIO.parse("data/input_sequences_rnd_9000.fasta", "fasta"))
    seq = str(rec.seq).upper()
    root_id = rec.id

    # 3) build the hierarchy:
    asm = HierarchicalAssembler(levels)
    root = asm.build(seq, root_id)

    # 4) export everything via export_utils:
    export_fragments_to_fasta(root, Path("data/tree_fragments.fasta"), root_id)
    export_tree_to_json(root, Path("data/tree_fragments.json"), root_id)
    export_tree_to_html(root, Path("data/tree_fragments.html"), root_id)
    export_clusters_and_fragments_html(
        root,
        Path("data/tree_clusters_fragments.html"),
        root_id
    )
