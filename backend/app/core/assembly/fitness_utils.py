# File: backend/app/core/assembly/fitness_utils.py
# Version: v0.1.3

"""
Fitness evaluation utilities for the GAOverlapSelector.

Provides:
  - Continuous orthogonality metrics (average and worst‐case normalized edit distances)
  - Quantitative penalties for 3′‐end mis‐annealing events
  - Quantitative penalties for k‐mer mis‐annealing hits
  - LRU‐caching of expensive operations (reverse complements, edit distances, k‐mer hashing)

These functions are designed to be called repeatedly in tight loops, so caching
and bit‐level optimizations (2‐bit encoding) are employed.
"""

import itertools
from functools import lru_cache
from typing import List, FrozenSet, Set

# === Configuration constants ===

# Minimum normalized edit distance threshold (unused in continuous scoring, but kept for reference)
EDIT_DISTANCE_MIN = 0.5

# Parameters for 3′‐end mis‐annealing detection
THREE_PRIME_LEN = 8
THREE_PRIME_MAX_MISMATCHES = 1

# k‐mer size for global mis‐annealing detection
KMER_SIZE = 10

# Two‐bit encoding for bases A, C, G, T
BASE_ENCODING = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}


@lru_cache(maxsize=None)
def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of `seq`.
    Cached to avoid recomputing for repeated sequences.
    """
    # Translate A<->T, C<->G then reverse
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@lru_cache(maxsize=None)
def encode_dna_2bit(seq: str) -> int:
    """
    Encode `seq` into a compact 2‐bit integer representation.
    A=00, C=01, G=10, T=11.
    Cached for reuse in k‐mer hashing and Hamming distance.
    """
    val = 0
    for base in seq:
        val = (val << 2) | BASE_ENCODING[base]
    return val


@lru_cache(maxsize=None)
def levenshtein_distance(a: str, b: str) -> int:
    """
    Compute the Levenshtein (edit) distance between strings a and b.
    Uses dynamic programming, with caching to reuse results.
    """
    # Ensure a is the longer
    if len(a) < len(b):
        return levenshtein_distance(b, a)
    # If the shorter string is empty, distance is length of the longer
    if not b:
        return len(a)

    # Initialize DP row
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
    """
    Compute the normalized edit distance: dist(a,b) / max(len(a), len(b)).
    Returns a value in [0.0, 1.0], where 1.0 means completely dissimilar.
    """
    # Special case: both empty
    if not a and not b:
        return 0.0
    dist = levenshtein_distance(a, b)
    return dist / max(len(a), len(b))


def average_normalized_distance(overlaps: List[str]) -> float:
    """
    Compute the average normalized edit distance over all pairwise combinations
    of the list `overlaps`.
    Returns 1.0 if fewer than 2 overlaps (maximal orthogonality).
    """
    pairs = list(itertools.combinations(overlaps, 2))
    if not pairs:
        return 1.0
    total = 0.0
    for a, b in pairs:
        total += normalized_edit_distance(a, b)
    return total / len(pairs)


def min_normalized_distance(overlaps: List[str]) -> float:
    """
    Compute the minimal (worst-case) normalized edit distance among all pairs
    of `overlaps`. Returns 1.0 if fewer than 2 overlaps.
    """
    pairs = list(itertools.combinations(overlaps, 2))
    if not pairs:
        return 1.0
    # Return the smallest normalized distance
    return min(normalized_edit_distance(a, b) for a, b in pairs)


@lru_cache(maxsize=None)
def extract_kmer_hashes(sequence: str, k: int) -> FrozenSet[int]:
    """
    Return the set of rolling k‐mer 2‐bit hashes for `sequence`.
    Uses a sliding‐window approach with bitmask to drop the oldest base.
    Cached to avoid re‐hashing the same substrings repeatedly.
    """
    if len(sequence) < k:
        return frozenset()
    hashes: Set[int] = set()

    # Encode the first k-mer
    curr_hash = encode_dna_2bit(sequence[:k])
    hashes.add(curr_hash)

    # Precompute mask to keep only 2*k bits
    mask = (1 << (2 * k)) - 1

    # Slide window by one base at a time
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
    """
    Count the number of potential mis‐annealing events where the reverse‐complement
    of the overlap's 3′ end (length three_prime_len) can bind within any of the
    sequences in `oligos`, with Hamming distance <= max_mismatches.
    Each sliding‐window hit counts as one event.
    """
    # Compute the reverse complement of the 3' tail
    rc_tail = reverse_complement(overlap[-three_prime_len:])
    tail_hash = encode_dna_2bit(rc_tail)
    mask      = (1 << (2 * three_prime_len)) - 1

    events = 0

    for oligo in oligos:
        if len(oligo) < three_prime_len:
            continue
        # Hash the first window
        window_hash = encode_dna_2bit(oligo[:three_prime_len])
        # Check initial window
        if hamming_distance_bits(window_hash, tail_hash, three_prime_len) <= max_mismatches:
            events += 1
        # Slide along the oligo
        for base in oligo[three_prime_len:]:
            window_hash = ((window_hash << 2) | BASE_ENCODING[base]) & mask
            if hamming_distance_bits(window_hash, tail_hash, three_prime_len) <= max_mismatches:
                events += 1

    return events


def count_rc_kmer_hits(
    overlap: str,
    oligos: List[str],
    k: int = KMER_SIZE
) -> int:
    """
    Count the total number of shared k‐mers between the reverse‐complement of `overlap`
    and all sequences in `oligos`. Each matching k‐mer counts once per occurrence.
    """
    rc_overlap = reverse_complement(overlap)
    ov_kmers   = extract_kmer_hashes(rc_overlap, k)
    total_hits = 0

    for oligo in oligos:
        # Intersection of k-mer sets yields shared k-mers
        total_hits += len(ov_kmers & extract_kmer_hashes(oligo, k))

    return total_hits


def hamming_distance_bits(x: int, y: int, length: int) -> int:
    """
    Compute the Hamming distance between two 2‐bit encoded sequences of `length`.
    Each pair of bits corresponds to one base; count mismatches by testing (xor & 0b11).
    """
    xor = x ^ y
    count = 0
    for _ in range(length):
        if (xor & 0b11) != 0:
            count += 1
        xor >>= 2
    return count
