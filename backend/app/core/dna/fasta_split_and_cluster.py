# File: backend/app/core/dna/fasta_split_and_cluster.py
# Version: v0.8.0

from typing import List, Tuple, Dict
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import json

from backend.app.core.assembly.ga_overlap_selector import GAOverlapSelector
from backend.app.core.visualization.oligo_html_visualizer import export_html_report, reverse_complement

# Parameters
OLIGO_LEN_MIN = 40
OLIGO_LEN_MAX = 60
NUM_OLIGO_MIN = 3
NUM_OLIGO_MAX = 6


class FastaParser:
    @staticmethod
    def parse_fasta(fpath: Path) -> Dict[str,str]:
        recs = SeqIO.parse(fpath, "fasta")
        return {r.id: str(r.seq).upper() for r in recs}


class OligoSplitter:
    def __init__(self, mn: int, mx: int):
        self.min_len = mn
        self.max_len = mx

    def split_sequence(self, seq: str) -> Tuple[List[str], List[Tuple[int,int]]]:
        bodies, poses, i = [], [], 0
        N = len(seq)
        while i < N:
            L = min(self.max_len, N - i)
            if L < self.min_len:
                break
            bodies.append(seq[i:i+L])
            poses.append((i, i+L))
            i += L
        return bodies, poses

    def assign_overlaps(
        self, full_seq: str,
        bodies: List[str],
        poses: List[Tuple[int,int]]
    ) -> Tuple[List[Tuple[str,str,int,int]], List[float]]:
        """
        Returns full-length oligos (seq, strand, start, end) and GA log.
        """
        n = len(bodies)
        if n < 2:
            s,e = poses[0]
            return ([(bodies[0], "+", s, e)], [])

        ga = GAOverlapSelector(full_sequence=full_seq,
                               oligo_positions=poses,
                               oligo_seqs=bodies)
        best, log = ga.evolve()
        ovs = best.overlaps  # List of (ov_seq, ov_start, ov_end)

        full_oligos: List[Tuple[str,str,int,int]] = []
        for idx, ((s0,e0), body) in enumerate(zip(poses, bodies)):
            prev_ov = ovs[idx-1] if idx>0 else None
            next_ov = ovs[idx]   if idx< n-1 else None

            start = prev_ov[1] if prev_ov else s0
            end   = next_ov[2] if next_ov else e0

            # build full sequence on sense
            frag = full_seq[start:end]
            strand = "+" if (idx % 2 == 0) else "-"
            seq = frag if strand=="+" else reverse_complement(frag)
            full_oligos.append((seq, strand, start, end))

        return full_oligos, log


def cluster_oligos(
    infos: List[Tuple[str,str,int,int]]
) -> List[List[Tuple[str,str,int,int]]]:
    clusters, i = [], 0
    M = len(infos)
    while i < M:
        sz = min(NUM_OLIGO_MAX, M - i)
        if sz < NUM_OLIGO_MIN:
            break
        clusters.append(infos[i:i+sz])
        i += sz
    return clusters


def export_oligos_to_fasta(
    clusters: Dict[str,List[List[Tuple[str,str,int,int]]]],
    outpath: Path
):
    recs = []
    for sid, cl_list in clusters.items():
        for ci, cl in enumerate(cl_list,1):
            for oi, (seq, strand, s, e) in enumerate(cl,1):
                rid = f"{sid}_c{ci}_o{oi}_{'sense' if strand=='+' else 'anti'}"
                recs.append(SeqRecord(Seq(seq), id=rid, description=""))
    SeqIO.write(recs, outpath, "fasta")


def export_oligos_to_json(
    clusters: Dict[str,List[List[Tuple[str,str,int,int]]]],
    outpath: Path
):
    out: Dict[str, List] = {}
    for sid, cl_list in clusters.items():
        out[sid] = []
        for ci, cl in enumerate(cl_list,1):
            info_list = []
            for oi, (seq, strand, s, e) in enumerate(cl,1):
                info_list.append({
                    "oligo_id": f"{sid}_c{ci}_o{oi}",
                    "sequence": seq,
                    "strand": strand,
                    "start": s,
                    "end": e,
                    "cluster_index": ci,
                    "oligo_index": oi
                })
            out[sid].append(info_list)
    with open(outpath, "w") as f:
        json.dump(out, f, indent=2)


def process_fasta_to_clusters(
    fasta_path: Path
) -> Tuple[
    Dict[str,str],
    Dict[str,List[List[Tuple[str,str,int,int]]]],
    Dict[str,List[float]]
]:
    full_seqs = FastaParser.parse_fasta(fasta_path)
    splitter  = OligoSplitter(OLIGO_LEN_MIN, OLIGO_LEN_MAX)

    cluster_map: Dict[str, List[List[Tuple[str,str,int,int]]]] = {}
    ga_logs:    Dict[str, List[float]] = {}

    for sid, seq in full_seqs.items():
        bodies, poses = splitter.split_sequence(seq)
        infos, log   = splitter.assign_overlaps(seq, bodies, poses)
        cluster_map[sid] = cluster_oligos(infos)
        ga_logs[sid]    = log

    return full_seqs, cluster_map, ga_logs


if __name__ == "__main__":
    fpath = Path("data/input_sequences.fasta")
    full, clusters, logs = process_fasta_to_clusters(fpath)
    export_oligos_to_fasta(clusters, Path("data/exported_oligos.fasta"))
    export_oligos_to_json(clusters, Path("data/exported_oligos.json"))
    export_html_report(
        Path("data/oligo_visualization.html"),
        clusters_by_sequence=clusters,
        ga_logs=logs,
        full_sequences=full
    )
