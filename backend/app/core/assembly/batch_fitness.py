# File: backend/app/core/assembly/batch_fitness.py
# Version: v0.1.0
"""
Batch & parallel helpers on top of fitness_utils.

Goals:
- Make GA evaluations fast and simple to call.
- Use threads (works well with edlib/primer3 C code that releases the GIL).
- Keep deterministic behavior (fixed seed is up to the caller).
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from math import comb
from typing import Iterable, List, Sequence, Tuple

from backend.app.core.assembly import fitness_utils as fu


@dataclass(slots=True)
class BatchOptions:
    """Controls parallelism and chunking."""
    workers: int = 0           # 0/None → auto = min(32, os.cpu_count() or 1)
    chunk_size: int = 128      # per-task inner size for small overhead
    tm_method: str = "PRIMER3" # default to fast primer3 if available


def _auto_workers(workers: int | None) -> int:
    if workers and workers > 0:
        return workers
    try:
        import os
        return max(1, min(32, os.cpu_count() or 1))
    except Exception:
        return 4


def _chunks(seq: Sequence, n: int) -> Iterable[Sequence]:
    if n <= 0:
        n = len(seq)
    for i in range(0, len(seq), n):
        yield seq[i : i + n]


# ------------------------
# Edit distance aggregates
# ------------------------

def pairwise_norm_edit_stats(
    overlaps: Sequence[str],
    *,
    options: BatchOptions | None = None,
) -> Tuple[float, float]:
    """
    Compute (avg_norm_edit_distance, min_norm_edit_distance) across all pairs.

    Parallelizes by splitting the combinations across workers.
    """
    opts = options or BatchOptions()
    m = len(overlaps)
    if m < 2:
        return 1.0, 1.0

    # Pre-materialize all pair indices (works well for m<=~1000). If m can be much larger,
    # switch to a streaming schedule. For CornStructor GA, m is typically small.
    pairs: List[Tuple[int, int]] = []
    pairs_extend = pairs.extend
    for i in range(m - 1):
        row = [(i, j) for j in range(i + 1, m)]
        pairs_extend(row)

    total_pairs = len(pairs)
    if total_pairs != comb(m, 2):
        # Guardrail; should never happen
        raise RuntimeError("pair enumeration mismatch")

    def eval_chunk(chunk: Sequence[Tuple[int, int]]) -> Tuple[float, float]:
        s = 0.0
        local_min = 1.0
        ov = overlaps
        for i, j in chunk:
            v = fu.normalized_edit_distance(ov[i], ov[j])
            s += v
            if v < local_min:
                local_min = v
        return s, local_min

    # Single-thread fast path
    if _auto_workers(opts.workers) == 1 or total_pairs <= opts.chunk_size:
        s, mm = eval_chunk(pairs)
        return s / total_pairs, mm

    # Parallel path
    agg_sum = 0.0
    global_min = 1.0
    with ThreadPoolExecutor(max_workers=_auto_workers(opts.workers)) as ex:
        futures = [ex.submit(eval_chunk, chunk) for chunk in _chunks(pairs, opts.chunk_size)]
        for f in as_completed(futures):
            s, mm = f.result()
            agg_sum += s
            if mm < global_min:
                global_min = mm
    return agg_sum / total_pairs, global_min


# -------------
# k-mer metrics
# -------------

def rc_kmer_hits_total(
    overlap: str,
    oligos: Sequence[str],
    k: int,
    *,
    options: BatchOptions | None = None,
) -> int:
    """
    Parallel total of shared k-mers between RC(overlap) and each oligo.
    """
    opts = options or BatchOptions()

    def eval_chunk(chunk: Sequence[str]) -> int:
        return sum(fu.count_rc_kmer_hits(overlap, [s], k=k) for s in chunk)

    if _auto_workers(opts.workers) == 1 or len(oligos) <= opts.chunk_size:
        return eval_chunk(list(oligos))

    total = 0
    with ThreadPoolExecutor(max_workers=_auto_workers(opts.workers)) as ex:
        futures = [ex.submit(eval_chunk, list(chunk)) for chunk in _chunks(oligos, opts.chunk_size)]
        for f in as_completed(futures):
            total += f.result()
    return total


# ------------
# Tm (batched)
# ------------

def compute_tm_batch(
    seqs: Sequence[str],
    *,
    method: str | None = None,
    options: BatchOptions | None = None,
) -> List[float]:
    """
    Compute Tm for many sequences, optionally in parallel.
    """
    opts = options or BatchOptions()
    meth = (method or opts.tm_method)

    def eval_chunk(chunk: Sequence[str]) -> List[float]:
        return [fu.compute_tm(s, method=meth) for s in chunk]

    if _auto_workers(opts.workers) == 1 or len(seqs) <= opts.chunk_size:
        return eval_chunk(list(seqs))

    results: List[float] = [0.0] * len(seqs)
    # Assign contiguous slices to avoid extra ordering logic
    with ThreadPoolExecutor(max_workers=_auto_workers(opts.workers)) as ex:
        futs = []
        start = 0
        for chunk in _chunks(seqs, opts.chunk_size):
            idx0 = start
            futs.append((idx0, ex.submit(eval_chunk, list(chunk))))
            start += len(chunk)
        for idx0, f in futs:
            vals = f.result()
            results[idx0: idx0 + len(vals)] = vals
    return results
