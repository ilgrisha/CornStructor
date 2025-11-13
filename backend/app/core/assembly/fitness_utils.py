# File: backend/app/core/assembly/fitness_utils.py
# Version: v0.4.2

"""
Fitness evaluation utilities for the GAOverlapSelector — C-accelerated.

Changes in v0.4.2:
- Add compatibility aliases: `gc_content_percent` -> `gc_fraction`,
  and `tm_of_sequence` -> `compute_tm`. This lets older/newer modules
  import either name without breaking.

Changes in v0.4.1:
- Switch from deprecated primer3.calcTm → primer3.calc_tm to silence warnings.
- Keep graceful fallbacks to Biopython NN/Wallace when primer3 is unavailable.
"""

from __future__ import annotations

import itertools
from functools import lru_cache
from typing import Iterable, FrozenSet, List, Set

# --- Optional accelerators ----------------------------------------------------

_EDLIB_AVAILABLE = False
_RAPIDFUZZ_AVAILABLE = False
_PRIMER3_AVAILABLE = False

try:
    import edlib  # type: ignore
    _EDLIB_AVAILABLE = True
except Exception:
    _EDLIB_AVAILABLE = False

if not _EDLIB_AVAILABLE:
    try:
        from rapidfuzz.distance import Levenshtein as _rf_lev  # type: ignore
        _RAPIDFUZZ_AVAILABLE = True
    except Exception:
        _RAPIDFUZZ_AVAILABLE = False

try:
    import primer3  # type: ignore
    _PRIMER3_AVAILABLE = True
except Exception:
    _PRIMER3_AVAILABLE = False

# Biopython fallback for Tm
try:
    from Bio.SeqUtils import MeltingTemp as _mt  # type: ignore
except Exception:  # pragma: no cover
    _mt = None  # type: ignore

# --- Configuration constants --------------------------------------------------

THREE_PRIME_LEN = 8
THREE_PRIME_MAX_MISMATCHES = 1
KMER_SIZE = 10

BASE_ENCODING = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}


# --- Small, hot helpers (cached) ---------------------------------------------

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


# --- Edit distance (C-accelerated with graceful fallbacks) -------------------

@lru_cache(maxsize=None)
def _lev_distance(a: str, b: str) -> int:
    """
    Edit (Levenshtein) distance using fastest available backend.

    Priority:
      1) edlib (task='distance')
      2) rapidfuzz.distance.Levenshtein.distance
      3) pure-Python fallback
    """
    if a is b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)

    if _EDLIB_AVAILABLE:
        return int(edlib.align(a, b, task="distance")["editDistance"])

    if _RAPIDFUZZ_AVAILABLE:
        return int(_rf_lev.distance(a, b))

    # Pure-Python fallback
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a):
        curr = [i + 1]
        for j, cb in enumerate(b):
            ins = prev[j + 1] + 1
            dele = curr[j] + 1
            sub = prev[j] + (ca != cb)
            curr.append(min(ins, dele, sub))
        prev = curr
    return prev[-1]


def normalized_edit_distance(a: str, b: str) -> float:
    """Normalized edit distance ∈ [0,1]."""
    if not a and not b:
        return 0.0
    dist = _lev_distance(a, b)
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


# --- Rolling k-mers (bit-packed) ---------------------------------------------

@lru_cache(maxsize=None)
def extract_kmer_hashes(sequence: str, k: int) -> FrozenSet[int]:
    """Set of rolling k-mer 2-bit hashes for `sequence` (cached)."""
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


# --- 3′-end mis-annealing checks --------------------------------------------

def hamming_distance_bits(x: int, y: int, length: int) -> int:
    """Hamming distance between 2-bit encoded sequences of given length."""
    xor = x ^ y
    count = 0
    for _ in range(length):
        if (xor & 0b11) != 0:
            count += 1
        xor >>= 2
    return count


def count_3prime_mismatch_events(
    overlap: str,
    oligos: List[str],
    three_prime_len: int = THREE_PRIME_LEN,
    max_mismatches: int = THREE_PRIME_MAX_MISMATCHES
) -> int:
    """Count reverse-complement 3′ tail binding events across oligos."""
    if three_prime_len <= 0 or not overlap:
        return 0
    rc_tail = reverse_complement(overlap[-three_prime_len:])
    tail_hash = encode_dna_2bit(rc_tail)
    mask = (1 << (2 * three_prime_len)) - 1
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
    """Total shared k-mers between RC(overlap) and all oligos."""
    if not overlap or k <= 0:
        return 0
    rc_overlap = reverse_complement(overlap)
    ov_kmers = extract_kmer_hashes(rc_overlap, k)
    total_hits = 0
    for oligo in oligos:
        total_hits += len(ov_kmers & extract_kmer_hashes(oligo, k))
    return total_hits


# --- Biophysical helpers -----------------------------------------------------

def gc_fraction(seq: str) -> float:
    """GC percent ∈ [0,100]."""
    if not seq:
        return 0.0
    s = seq.upper()
    g = s.count('G')
    c = s.count('C')
    return (g + c) / len(s) * 100.0


def max_same_base_run(seq: str) -> int:
    """Maximum homopolymer run length in the sequence."""
    if not seq:
        return 0
    best = 1
    curr = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            curr += 1
            if curr > best:
                best = curr
        else:
            curr = 1
    return best


def contains_any_motif(seq: str, motifs: Iterable[str]) -> bool:
    """True if sequence contains any exact motif (case-insensitive)."""
    if not seq:
        return False
    s = seq.upper()
    return any(m and (m.upper() in s) for m in motifs)


def count_motif_occurrences(seq: str, motif: str) -> int:
    """Count overlapping exact motif occurrences (case-insensitive)."""
    if not seq or not motif:
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


# --- Melting temperature (primer3 preferred) ---------------------------------

def compute_tm(seq: str, method: str = "PRIMER3", **params) -> float:
    """
    Compute melting temperature (°C).

    method:
      - "PRIMER3" (default): via `primer3.calc_tm`
          mv_conc  = K + Na (mM)
          dv_conc  = Mg (mM)
          dntp_conc= dNTPs (mM)
          dna_conc = max(dnac1, dnac2) (nM)
      - "NN": Biopython nearest-neighbor
      - "Wallace": Biopython Wallace rule

    Falls back automatically if the requested method is unavailable.
    """
    seq = (seq or "").upper()

    # Q5-like defaults
    Na_mM = float(params.get("Na", 0.0))
    K_mM = float(params.get("K", 50.0))
    Tris_mM = float(params.get("Tris", 10.0))
    Mg_mM = float(params.get("Mg", 2.0))
    dNTPs_mM = float(params.get("dNTPs", 0.8))
    dnac1_nM = float(params.get("dnac1", 500.0))
    dnac2_nM = float(params.get("dnac2", 500.0))
    saltcorr = int(params.get("saltcorr", 7))

    method_u = (method or "").upper()

    if method_u in ("PRIMER3", "P3", "PR3"):
        if _PRIMER3_AVAILABLE:
            mv_conc = K_mM + Na_mM
            dv_conc = Mg_mM
            dna_conc = max(dnac1_nM, dnac2_nM)
            # Use the non-deprecated API
            return float(
                primer3.calc_tm(
                    seq,
                    mv_conc=mv_conc,
                    dv_conc=dv_conc,
                    dntp_conc=dNTPs_mM,
                    dna_conc=dna_conc,
                )
            )
        # fall through to Biopython if primer3 missing

    if method_u == "WALLACE":
        if _mt is None:
            raise RuntimeError("Biopython MeltingTemp not available.")
        return float(_mt.Tm_Wallace(seq))

    # NN (nearest neighbor) via Biopython
    if _mt is None:
        raise RuntimeError("Biopython MeltingTemp not available and primer3 disabled.")
    return float(
        _mt.Tm_NN(
            seq,
            Na=Na_mM,
            K=K_mM,
            Tris=Tris_mM,
            Mg=Mg_mM,
            dNTPs=dNTPs_mM,
            dnac1=dnac1_nM,
            dnac2=dnac2_nM,
            saltcorr=saltcorr,
        )
    )


# --- Compatibility aliases ----------------------------------------------------
def gc_content_percent(seq: str) -> float:
    """
    Back-compat alias used by newer planner/GA modules.
    Equivalent to gc_fraction(seq).
    """
    return gc_fraction(seq)


def tm_of_sequence(seq: str, method: str = "PRIMER3", **params) -> float:
    """
    Back-compat alias used by newer planner/GA modules.
    Equivalent to compute_tm(seq, method=..., **params).
    """
    return compute_tm(seq, method=method, **params)
