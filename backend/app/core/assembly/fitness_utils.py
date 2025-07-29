# File: backend/app/core/assembly/fitness_utils.py
# Version: v0.1.0

"""
Fitness evaluation utilities for genetic algorithm.
Includes orthogonality, 3' misannealing, and global k-mer misannealing checks.
"""

import itertools
from typing import List

EDIT_DISTANCE_MIN = 0.5
THREE_PRIME_LEN = 8
THREE_PRIME_MAX_MISMATCHES = 1
KMER_SIZE = 10
MAX_KMER_HITS = 3

BASE_ENCODING = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def encode_dna_2bit(seq: str) -> int:
    val = 0
    for base in seq:
        val = (val << 2) | BASE_ENCODING[base]
    return val


def levenshtein_distance(a: str, b: str) -> int:
    if len(a) < len(b):
        return levenshtein_distance(b, a)
    if len(b) == 0:
        return len(a)
    prev_row = list(range(len(b) + 1))
    for i, c1 in enumerate(a):
        curr_row = [i + 1]
        for j, c2 in enumerate(b):
            insertions = prev_row[j + 1] + 1
            deletions = curr_row[j] + 1
            substitutions = prev_row[j] + (c1 != c2)
            curr_row.append(min(insertions, deletions, substitutions))
        prev_row = curr_row
    return prev_row[-1]


def is_orthogonal(overlaps: List[str], normalized_min: float = EDIT_DISTANCE_MIN) -> bool:
    for a, b in itertools.combinations(overlaps, 2):
        max_len = max(len(a), len(b))
        dist = levenshtein_distance(a, b)
        norm = dist / max_len
        if norm < normalized_min:
            return False
    return True


def has_bad_3prime_binding_fast(overlap: str, oligos: List[str],
                                three_prime_len: int = THREE_PRIME_LEN,
                                max_mismatches: int = THREE_PRIME_MAX_MISMATCHES) -> bool:
    rc_3prime = reverse_complement(overlap[-three_prime_len:])
    encoded_query = encode_dna_2bit(rc_3prime)
    mask = (1 << (2 * three_prime_len)) - 1

    for oligo in oligos:
        if len(oligo) < three_prime_len:
            continue
        encoded_oligo = encode_dna_2bit(oligo[:three_prime_len])
        if hamming_distance_bits(encoded_query, encoded_oligo, three_prime_len) <= max_mismatches:
            return True
        for i in range(three_prime_len, len(oligo)):
            encoded_oligo = ((encoded_oligo << 2) | BASE_ENCODING[oligo[i]]) & mask
            if hamming_distance_bits(encoded_query, encoded_oligo, three_prime_len) <= max_mismatches:
                return True
    return False


def extract_kmer_hashes(sequence: str, k: int) -> set:
    if len(sequence) < k:
        return set()
    hashes = set()
    curr_hash = encode_dna_2bit(sequence[:k])
    hashes.add(curr_hash)
    mask = (1 << (2 * k)) - 1
    for i in range(k, len(sequence)):
        curr_hash = ((curr_hash << 2) | BASE_ENCODING[sequence[i]]) & mask
        hashes.add(curr_hash)
    return hashes


def rc_kmer_misannealing(overlap: str, oligos: List[str],
                          k: int = KMER_SIZE, max_hits: int = MAX_KMER_HITS) -> bool:
    rc_overlap = reverse_complement(overlap)
    overlap_kmers = extract_kmer_hashes(rc_overlap, k)
    hit_count = 0
    for oligo in oligos:
        oligo_kmers = extract_kmer_hashes(oligo, k)
        matches = overlap_kmers & oligo_kmers
        hit_count += len(matches)
        if hit_count >= max_hits:
            return True
    return False


def hamming_distance_bits(x: int, y: int, length: int) -> int:
    xor = x ^ y
    count = 0
    for _ in range(length):
        if (xor & 0b11) != 0:
            count += 1
        xor >>= 2
    return count
