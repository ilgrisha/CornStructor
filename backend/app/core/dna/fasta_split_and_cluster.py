# File: backend/app/core/dna/fasta_split_and_cluster.py
# Version: v0.8.4

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

        Args:
            fpath: Path to the FASTA file.

        Returns:
            Dict of {sequence_id: sequence_str}.
        """
        records = SeqIO.parse(fpath, "fasta")
        return {rec.id: str(rec.seq).upper() for rec in records}


class OligoSplitter:
    """
    Splits a long DNA sequence into oligo 'bodies' and uses a GA to choose overlaps.
    Reconstructs full-length oligo sequences (with overlaps) and assigns strand orientation.
    """
    def __init__(self, min_len: int, max_len: int):
        """
        Args:
            min_len: Minimum length of each oligo body (without overlaps).
            max_len: Maximum length of each oligo body.
        """
        self.min_len = min_len
        self.max_len = max_len

    def split_sequence(self, seq: str) -> Tuple[List[str], List[Tuple[int, int]]]:
        """
        Greedy windowing: split `seq` into non-overlapping oligo bodies.

        Args:
            seq: Full DNA sequence to split.

        Returns:
            bodies: List of oligo substrings.
            poses:  List of (start, end) positions of each body in `seq`.
        """
        bodies: List[str] = []
        poses: List[Tuple[int, int]] = []
        i = 0
        N = len(seq)

        # Slide window until end of sequence
        while i < N:
            L = min(self.max_len, N - i)
            if L < self.min_len:
                # Remaining tail too short to form an oligo
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
        """
        Use GAOverlapSelector to pick optimal overlaps at each junction, then build
        full-length oligos (including upstream and downstream overlaps),
        alternating strand orientation (+ for sense, - for antisense).

        Args:
            full_seq: The original full-length DNA sequence.
            bodies:   List of oligo bodies (non-overlapping substrings).
            poses:    Corresponding (start,end) positions for each body.

        Returns:
            full_oligos: List of tuples
              (full_seq, strand, start, end, overlap_prev, overlap_next)
            ga_log:      GA fitness log (best score per generation).
        """
        n = len(bodies)
        # If only one body, just return it (no overlaps needed)
        if n < 2:
            s0, e0 = poses[0]
            return ([(bodies[0], "+", s0, e0, "", "")], [])

        # Initialize and run the GA to choose overlaps
        ga = GAOverlapSelector(
            full_sequence=full_seq,
            oligo_positions=poses,
            oligo_seqs=bodies
        )
        best, ga_log = ga.evolve()
        # best.overlaps is List[ (ov_seq, ov_start, ov_end) ]
        ovs = best.overlaps

        full_oligos: List[Tuple[str, str, int, int, str, str]] = []
        for idx, ((s0, e0), body) in enumerate(zip(poses, bodies)):
            # Determine upstream/downstream overlap metadata
            prev_ov = ovs[idx - 1] if idx > 0 else None
            next_ov = ovs[idx]     if idx < n - 1 else None

            # Compute full-oligo genomic bounds
            start = prev_ov[1] if prev_ov else s0
            end   = next_ov[2] if next_ov else e0
            frag  = full_seq[start:end]

            # Alternate strand: even idx → sense (+), odd idx → antisense (-)
            strand   = "+" if (idx % 2 == 0) else "-"
            seq_full = frag if strand == "+" else reverse_complement(frag)

            # Extract raw overlap sequences
            prev_seq = prev_ov[0] if prev_ov else ""
            next_seq = next_ov[0]  if next_ov else ""

            full_oligos.append((seq_full, strand, start, end, prev_seq, next_seq))

        return full_oligos, ga_log


def cluster_oligos(
    infos: List[Tuple[str, str, int, int, str, str]]
) -> List[List[Tuple[str, str, int, int, str, str]]]:
    """
    Group a flat list of oligo infos into clusters of size between NUM_OLIGO_MIN and NUM_OLIGO_MAX.

    Args:
        infos: Flat list of (seq, strand, start, end, prev_ov, next_ov).

    Returns:
        List of clusters; each cluster is a list of the above tuples.
    """
    clusters: List[List[Tuple[str, str, int, int, str, str]]] = []
    i = 0
    M = len(infos)

    # Greedily take up to NUM_OLIGO_MAX per cluster
    while i < M:
        sz = min(NUM_OLIGO_MAX, M - i)
        if sz < NUM_OLIGO_MIN:
            # If remaining bodies < minimum cluster size, stop
            break
        clusters.append(infos[i : i + sz])
        i += sz

    return clusters


def export_oligos_to_fasta(
    clusters: Dict[str, List[List[Tuple[str, str, int, int, str, str]]]],
    outpath: Path
):
    """
    Write all full-length oligos (including overlaps) to a multi-FASTA file.
    ID format: <seq_id>_c<cluster>_o<oligo>_<sense|antisense>
    """
    records: List[SeqRecord] = []

    for sid, cl_list in clusters.items():
        for ci, cl in enumerate(cl_list, start=1):
            for oi, (seq, strand, s, e, prev_ov, next_ov) in enumerate(cl, start=1):
                label = "sense" if strand == "+" else "antisense"
                rid   = f"{sid}_c{ci}_o{oi}_{label}"
                records.append(SeqRecord(Seq(seq), id=rid, description=""))

    SeqIO.write(records, outpath, "fasta")


def export_oligos_to_json(
    clusters: Dict[str, List[List[Tuple[str, str, int, int, str, str]]]],
    outpath: Path
):
    """
    Export full-length oligo metadata to JSON, including:
      - oligo_id, sequence, strand, start, end
      - overlap_prev, overlap_next (reverse-complemented for antisense)
      - cluster & oligo indices
    """
    out: Dict[str, List] = {}

    for sid, cl_list in clusters.items():
        out[sid] = []
        for ci, cl in enumerate(cl_list, start=1):
            info_list = []
            for oi, (seq, strand, s, e, prev_ov, next_ov) in enumerate(cl, start=1):
                # Adjust overlap orientation if antisense
                if strand == "-":
                    prev_rc = reverse_complement(prev_ov) if prev_ov else ""
                    next_rc = reverse_complement(next_ov) if next_ov else ""
                else:
                    prev_rc = prev_ov
                    next_rc = next_ov

                info = {
                    "oligo_id":      f"{sid}_c{ci}_o{oi}",
                    "sequence":      seq,
                    "strand":        strand,
                    "start":         s,
                    "end":           e,
                    "overlap_prev":  prev_rc,
                    "overlap_next":  next_rc,
                    "cluster_index": ci,
                    "oligo_index":   oi
                }
                info_list.append(info)
            out[sid].append(info_list)

    # Write JSON with pretty indentation
    with open(outpath, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)


def process_fasta_to_clusters(
    fasta_path: Path
) -> Tuple[
    Dict[str, str],
    Dict[str, List[List[Tuple[str, str, int, int, str, str]]]],
    Dict[str, List[float]]
]:
    """
    High-level pipeline:
      1. Parse FASTA into sequences.
      2. Split each into oligo bodies.
      3. Run GA to assign overlaps and reconstruct full oligos.
      4. Cluster oligos.
      5. Return `full_sequences`, `cluster_map`, and GA `logs`.
    """
    # 1) Read input FASTA
    full_seqs = FastaParser.parse_fasta(fasta_path)

    # 2) Prepare the splitter
    splitter = OligoSplitter(OLIGO_LEN_MIN, OLIGO_LEN_MAX)

    cluster_map: Dict[str, List[List[Tuple[str, str, int, int, str, str]]]] = {}
    ga_logs:    Dict[str, List[float]]                                = {}

    # Process each sequence independently
    for sid, seq in full_seqs.items():
        bodies, poses       = splitter.split_sequence(seq)
        infos, log          = splitter.assign_overlaps(seq, bodies, poses)
        cluster_map[sid]    = cluster_oligos(infos)
        ga_logs[sid]        = log

    return full_seqs, cluster_map, ga_logs


# -----------------------------------------------------------------------------
# Example usage when run as a script
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    fasta_file = Path("data/input_sequences.fasta")
    full, clusters, logs = process_fasta_to_clusters(fasta_file)

    # Export designed oligos
    export_oligos_to_fasta(clusters, Path("data/exported_oligos.fasta"))
    export_oligos_to_json(clusters, Path("data/exported_oligos.json"))
    export_html_report(
        Path("data/oligo_visualization.html"),
        clusters_by_sequence=clusters,
        ga_logs=logs,
        full_sequences=full
    )
