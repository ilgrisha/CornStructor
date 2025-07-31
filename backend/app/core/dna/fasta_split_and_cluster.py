# File: backend/app/core/dna/fasta_split_and_cluster.py
# Version: v0.8.5

"""
Parses a FASTA file and splits each DNA sequence into clusters of overlapping oligos.
Uses a Genetic Algorithm (GAOverlapSelector) to select optimal overlap segments at each junction,
then reconstructs full-length oligo sequences (overlaps + body + overlaps), alternating strand orientation.
Exports designed oligos to:
  - FASTA (for synthesis)
  - JSON (with metadata including overlaps)
  - HTML (interactive visualization with aligned strands and GA progress)
"""

from typing import List, Tuple, Dict
from pathlib import Path
import json

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from backend.app.core.assembly.ga_overlap_selector import GAOverlapSelector
from backend.app.core.visualization.oligo_html_visualizer import (
    export_html_report,
    reverse_complement
)

# -----------------------------------------------------------------------------
# User-configurable parameters for oligo design
# -----------------------------------------------------------------------------
OLIGO_LEN_MIN = 40    # minimum oligo body length (bp)
OLIGO_LEN_MAX = 60    # maximum oligo body length (bp)
NUM_OLIGO_MIN = 3     # minimum oligos per cluster
NUM_OLIGO_MAX = 6     # maximum oligos per cluster


class FastaParser:
    """
    Utility class to parse FASTA files.
    """
    @staticmethod
    def parse_fasta(fpath: Path) -> Dict[str, str]:
        """
        Read a standard multiline FASTA file and return a dict mapping record IDs to uppercase sequences.
        """
        records = SeqIO.parse(fpath, "fasta")
        return {rec.id: str(rec.seq).upper() for rec in records}


class OligoSplitter:
    """
    Splits a long DNA sequence into oligo 'bodies' and uses a GA to choose overlaps.
    Reconstructs full-length oligo sequences (with overlaps) and assigns strand orientation.
    """
    def __init__(self, min_len: int, max_len: int):
        self.min_len = min_len
        self.max_len = max_len

    def split_sequence(self, seq: str) -> Tuple[List[str], List[Tuple[int, int]]]:
        bodies, poses, i = [], [], 0
        N = len(seq)
        while i < N:
            L = min(self.max_len, N - i)
            if L < self.min_len:
                break
            bodies.append(seq[i : i + L])
            poses.append((i, i + L))
            i += L
        return bodies, poses

    def assign_overlaps(
        self,
        full_seq: str,
        bodies: List[str],
        poses: List[Tuple[int, int]]
    ) -> Tuple[List[Tuple[str, str, int, int, str, str]], List[float]]:
        n = len(bodies)
        if n < 2:
            s0, e0 = poses[0]
            return ([(bodies[0], "+", s0, e0, "", "")], [])

        ga = GAOverlapSelector(
            full_sequence    = full_seq,
            oligo_positions  = poses,
            oligo_seqs       = bodies
        )
        best, ga_log = ga.evolve()
        ovs = best.overlaps

        full_oligos: List[Tuple[str, str, int, int, str, str]] = []
        for idx, ((s0, e0), body) in enumerate(zip(poses, bodies)):
            prev = ovs[idx-1] if idx>0      else None
            nxt  = ovs[idx]   if idx< n-1   else None

            start    = prev[1] if prev else s0
            end      = nxt[2]  if nxt  else e0
            frag     = full_seq[start:end]
            strand   = "+" if idx%2==0     else "-"
            seq_full = frag if strand=="+" else reverse_complement(frag)
            prev_seq = prev[0] if prev else ""
            next_seq = nxt[0]  if nxt  else ""

            full_oligos.append((seq_full, strand, start, end, prev_seq, next_seq))

        return full_oligos, ga_log


def cluster_oligos(
    infos: List[Tuple[str, str, int, int, str, str]]
) -> List[List[Tuple[str, str, int, int, str, str]]]:
    clusters, i, M = [], 0, len(infos)
    while i < M:
        sz = min(NUM_OLIGO_MAX, M - i)
        if sz < NUM_OLIGO_MIN:
            break
        clusters.append(infos[i : i + sz])
        i += sz
    return clusters


def export_oligos_to_fasta(
    clusters: Dict[str, List[List[Tuple[str, str, int, int, str, str]]]],
    outpath: Path
):
    from Bio.SeqRecord import SeqRecord
    records = []
    for sid, cl_list in clusters.items():
        for ci, cl in enumerate(cl_list, start=1):
            for oi, (seq, strand, s, e, prev, nxt) in enumerate(cl, start=1):
                label = "sense" if strand=="+" else "antisense"
                rid   = f"{sid}_c{ci}_o{oi}_{label}"
                records.append(SeqRecord(Seq(seq), id=rid, description=""))
    SeqIO.write(records, outpath, "fasta")


def export_oligos_to_json(
    clusters: Dict[str, List[List[Tuple[str, str, int, int, str, str]]]],
    outpath: Path
):
    out: Dict[str, List] = {}
    for sid, cl_list in clusters.items():
        out[sid] = []
        for ci, cl in enumerate(cl_list, start=1):
            lst = []
            for oi, (seq, strand, s, e, prev, nxt) in enumerate(cl, start=1):
                if strand=="-":
                    prev_rc = reverse_complement(prev) if prev else ""
                    nxt_rc  = reverse_complement(nxt)  if nxt  else ""
                else:
                    prev_rc, nxt_rc = prev, nxt
                lst.append({
                    "oligo_id":     f"{sid}_c{ci}_o{oi}",
                    "sequence":     seq,
                    "strand":       strand,
                    "start":        s,
                    "end":          e,
                    "overlap_prev": prev_rc,
                    "overlap_next": nxt_rc,
                    "cluster_index":ci,
                    "oligo_index":  oi
                })
            out[sid].append(lst)
    with open(outpath, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)


def process_fasta_to_clusters(
    fasta_path: Path
):
    full_seqs = FastaParser.parse_fasta(fasta_path)
    splitter  = OligoSplitter(OLIGO_LEN_MIN, OLIGO_LEN_MAX)
    cluster_map, ga_logs = {}, {}

    for sid, seq in full_seqs.items():
        bodies, poses       = splitter.split_sequence(seq)
        infos, log          = splitter.assign_overlaps(seq, bodies, poses)
        cluster_map[sid]    = cluster_oligos(infos)
        ga_logs[sid]        = log

    return full_seqs, cluster_map, ga_logs


if __name__ == "__main__":
    fasta_file = Path("data/input_sequences_rnd_700.fasta")
    full, clusters, logs = process_fasta_to_clusters(fasta_file)

    export_oligos_to_fasta(clusters, Path("data/exported_oligos.fasta"))
    export_oligos_to_json(clusters, Path("data/exported_oligos.json"))
    export_html_report(
        Path("data/oligo_visualization.html"),
        clusters_by_sequence=clusters,
        ga_logs=logs,
        full_sequences=full
    )
