# File: backend/app/core/assembly/fitness_utils.py
# Version: v0.3.1

"""
Fitness evaluation utilities for the GAOverlapSelector.

Provides:
  - Continuous orthogonality metrics (average and worst‐case normalized edit distances)
  - Quantitative penalties for 3′‐end mis‐annealing events
  - Quantitative penalties for k‐mer mis‐annealing hits
  - Biophysical helpers: GC fraction, homopolymer runs, motif checks and counts
  - Biopython-backed Tm computation (NN or Wallace)
  - LRU‐caching of expensive operations (reverse complements, edit distances, k‐mer hashing)

Defaults (v0.3.1):
  - Tm defaults approximate **Q5® 1× PCR** conditions:
      Na = 0.0 mM
      K  = 50.0 mM
      Tris = 10.0 mM
      Mg = 2.0 mM
      dNTPs = 0.8 mM   (200 µM each)
      Primer strand concentrations: dnac1 = dnac2 = 500 nM
      saltcorr = 7
"""

import itertools
from functools import lru_cache
from typing import List, FrozenSet, Set, Iterable, Optional

from Bio.SeqUtils import MeltingTemp as mt

# === Configuration constants ===

# Parameters for 3′‐end mis‐annealing detection
THREE_PRIME_LEN = 8
THREE_PRIME_MAX_MISMATCHES = 1

# k‐mer size for global mis‐annealing detection
KMER_SIZE = 10

# Two‐bit encoding for bases A, C, G, T
BASE_ENCODING = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}


@lru_cache(maxsize=None)
def reverse_complement(seq: str) -> str:
    """Reverse complement (cached)."""
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


@lru_cache(maxsize=None)
def encode_dna_2bit(seq: str) -> int:
    """2-bit encode sequence (cached)."""
    val = 0
    for base in seq:
        val = (val << 2) | BASE_ENCODING[base]
    return val


@lru_cache(maxsize=None)
def levenshtein_distance(a: str, b: str) -> int:
    """Levenshtein edit distance (cached)."""
    if len(a) < len(b):
        return levenshtein_distance(b, a)
    if not b:
        return len(a)
    prev_row = list(range(len(b) + 1))
    for i, ca in enumerate(a):
        curr_row = [i + 1]
        for j, cb in enumerate(b):
            ins  = prev_row[j + 1] + 1
            dele = curr_row[j]     + 1
            sub  = prev_row[j]     + (ca != cb)
            curr_row.append(min(ins, dele, sub))
        prev_row = curr_row
    return prev_row[-1]


def normalized_edit_distance(a: str, b: str) -> float:
    """Normalized edit distance ∈ [0,1]."""
    if not a and not b:
        return 0.0
    dist = levenshtein_distance(a, b)
    return dist / max(len(a), len(b))


def average_normalized_distance(overlaps: List[str]) -> float:
    """Average normalized edit distance across all pairs (1.0 if <2 overlaps)."""
    pairs = list(itertools.combinations(overlaps, 2))
    if not pairs:
        return 1.0
    total = 0.0
    for x, y in pairs:
        total += normalized_edit_distance(x, y)
    return total / len(pairs)


def min_normalized_distance(overlaps: List[str]) -> float:
    """Minimal pairwise normalized edit distance (1.0 if <2 overlaps)."""
    pairs = list(itertools.combinations(overlaps, 2))
    if not pairs:
        return 1.0
    return min(normalized_edit_distance(a, b) for a, b in pairs)


@lru_cache(maxsize=None)
def extract_kmer_hashes(sequence: str, k: int) -> FrozenSet[int]:
    """Set of rolling k‐mer 2‐bit hashes for `sequence` (cached)."""
    if len(sequence) < k:
        return frozenset()
    hashes: Set[int] = set()
    curr_hash = encode_dna_2bit(sequence[:k])
    hashes.add(curr_hash)
    mask = (1 << (2 * k)) - 1
    for base in sequence[k:]:
        curr_hash = ((curr_hash << 2) | BASE_ENCODING[base]) & mask
        hashes.add(curr_hash)
    return frozenset(hashes)


def count_3prime_mismatch_events(
    overlap: str,
    oligos: List[str],
    three_prime_len: int = THREE_PRIME_LEN,
    max_mismatches: int = THREE_PRIME_MAX_MISMATCHES
) -> int:
    """Count reverse‐complement 3′ tail binding events across oligos with Hamming ≤ threshold."""
    rc_tail = reverse_complement(overlap[-three_prime_len:])
    tail_hash = encode_dna_2bit(rc_tail)
    mask      = (1 << (2 * three_prime_len)) - 1
    events = 0
    for oligo in oligos:
        if len(oligo) < three_prime_len:
            continue
        window_hash = encode_dna_2bit(oligo[:three_prime_len])
        if hamming_distance_bits(window_hash, tail_hash, three_prime_len) <= max_mismatches:
            events += 1
        for base in oligo[three_prime_len:]:
            window_hash = ((window_hash << 2) | BASE_ENCODING[base]) & mask
            if hamming_distance_bits(window_hash, tail_hash, three_prime_len) <= max_mismatches:
                events += 1
    return events


def count_rc_kmer_hits(overlap: str, oligos: List[str], k: int = KMER_SIZE) -> int:
    """Total shared k‐mers between RC(overlap) and oligos."""
    rc_overlap = reverse_complement(overlap)
    ov_kmers   = extract_kmer_hashes(rc_overlap, k)
    total_hits = 0
    for oligo in oligos:
        total_hits += len(ov_kmers & extract_kmer_hashes(oligo, k))
    return total_hits


def hamming_distance_bits(x: int, y: int, length: int) -> int:
    """Hamming distance between 2‐bit encoded sequences of given length."""
    xor = x ^ y
    count = 0
    for _ in range(length):
        if (xor & 0b11) != 0:
            count += 1
        xor >>= 2
    return count


# ---------------------------
# New helpers for constraints
# ---------------------------

def gc_fraction(seq: str) -> float:
    """GC persent ∈ [0,100]."""
    if not seq:
        return 0.0
    g = seq.count('G')
    c = seq.count('C')
    return (g + c) / len(seq) * 100


def max_same_base_run(seq: str) -> int:
    """Maximum homopolymer run length in the sequence."""
    if not seq:
        return 0
    best = 1
    curr = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            curr += 1
            best = max(best, curr)
        else:
            curr = 1
    return best


def contains_any_motif(seq: str, motifs: Iterable[str]) -> bool:
    """True if sequence contains any exact motif."""
    s = seq.upper()
    return any(m and (m.upper() in s) for m in motifs)


def count_motif_occurrences(seq: str, motif: str) -> int:
    """Count overlapping exact motif occurrences."""
    if not motif:
        return 0
    s = seq.upper()
    m = motif.upper()
    count = start = 0
    while True:
        idx = s.find(m, start)
        if idx == -1:
            break
        count += 1
        start = idx + 1
    return count


def compute_tm(seq: str, method: str = "NN", **params) -> float:
    """
    Compute melting temperature using Biopython.

    Args:
        seq: DNA sequence (string)
        method:
          - "NN"      → mt.Tm_NN (nearest-neighbor)
          - "Wallace" → mt.Tm_Wallace
        params (for "NN"):
          Ion concentrations in **mM**: Na, K, Tris, Mg, dNTPs
          Strand concentrations in **nM**: dnac1, dnac2
          saltcorr: integer (recommended 7)

    Returns:
        Temperature in °C (float).

    Notes:
        This wrapper expects Q5-like defaults in mM/nM and forwards them directly
        to Biopython, which accepts these units for Tm_NN.
    """
    seq = seq.upper()
    if method.upper() == "WALLACE":
        return float(mt.Tm_Wallace(seq))

    # --- Q5 1× defaults (mM / nM) ---
    Na_mM   = float(params.get("Na", 0.0))
    K_mM    = float(params.get("K", 50.0))
    Tris_mM = float(params.get("Tris", 10.0))
    Mg_mM   = float(params.get("Mg", 2.0))
    dNTPs_mM= float(params.get("dNTPs", 0.8))
    dnac1_nM= float(params.get("dnac1", 500.0))
    dnac2_nM= float(params.get("dnac2", 500.0))
    saltcorr= int(params.get("saltcorr", 7))

    return float(mt.Tm_NN(
        seq,
        Na=Na_mM,
        K=K_mM,
        Tris=Tris_mM,
        Mg=Mg_mM,
        dNTPs=dNTPs_mM,
        dnac1=dnac1_nM,
        dnac2=dnac2_nM,
        saltcorr=saltcorr
    ))
