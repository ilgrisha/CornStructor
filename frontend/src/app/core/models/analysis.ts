// File: frontend/src/app/core/models/analysis.ts
// Version: v0.2.0

/**
 * Shared types for sequence analysis overlays and parameters.
 *
 * This version adds support for secondary-structure stems (ViennaRNA-backed),
 * keeping full backward compatibility with v0.1.0:
 *  - Extends FeatureKey with 'stems'
 *  - Extends FeatureToggles with `stems`
 *  - Extends AnalysisResult with `stems: Region[]`
 *  - Extends AnalysisParams with `stemMinLen` and `stemMergeMaxGap`
 *
 * Notes:
 *  - Indexing is 0-based, [start, end) half-open intervals.
 *  - Region remains generic (no `kind` field) to preserve existing renderers.
 */

export type FeatureKey = 'homopolymers' | 'gc' | 'entropy' | 'repeats' | 'stems';

export interface Region {
  /** 0-based inclusive start index */
  start: number;
  /** 0-based exclusive end index */
  end: number;
  /** Optional label for UI (e.g., motif, score, etc.) */
  label?: string;
  /** Arbitrary metadata bag for advanced overlays */
  extra?: Record<string, unknown>;
}

export interface AnalysisParams {
  // === Homopolymers ===
  /** Minimal homopolymer length to flag */
  homoMin: number;

  // === GC% windows ===
  /** Sliding window size for GC% computation */
  gcWindow: number;
  /** Lower GC% bound (0..100) to flag as out-of-range */
  gcMin: number;
  /** Upper GC% bound (0..100) to flag as out-of-range */
  gcMax: number;

  // === Low complexity via Shannon entropy ===
  /** Sliding window size for entropy calculation */
  entropyWindow: number;
  /** Entropy threshold in bits (0..2 on A/C/G/T) to flag as low complexity */
  entropyMin: number;

  // === Simple tandem repeats (non-homopolymer motifs) ===
  /** Minimum motif length (in bp) for tandem repeat scanning (≥2 to exclude homopolymers) */
  repeatKmin: number;
  /** Maximum motif length (in bp) for tandem repeat scanning */
  repeatKmax: number;
  /** Minimal number of consecutive motif repetitions to flag */
  repeatMinReps: number;

  // === NEW: Secondary structure stems (ViennaRNA-backed) ===
  /**
   * Minimal number of sequential nucleotides in a stem *after merging*
   * (i.e., the final linear region length threshold).
   */
  stemMinLen: number;
  /**
   * Merge adjacent stem regions if the number of unpaired nucleotides
   * between them is <= this value.
   */
  stemMergeMaxGap: number;
}

export interface FeatureToggles {
  homopolymers: boolean;
  gc: boolean;
  entropy: boolean;
  repeats: boolean;
  /** NEW: toggle visualization of stems */
  stems: boolean;
}

export interface AnalysisResult {
  /** Runs of identical bases length ≥ homoMin */
  homopolymers: Region[];
  /** Regions where GC% is outside [gcMin, gcMax] over window gcWindow */
  gcOutOfRange: Region[];
  /** Low-complexity bands where Shannon entropy < entropyMin over window entropyWindow */
  lowComplexity: Region[];
  /** Simple tandem repeats (motif length in [repeatKmin, repeatKmax], reps ≥ repeatMinReps) */
  repeats: Region[];
  /** NEW: Secondary-structure stems (merged linear regions of paired nts) */
  stems: Region[];
}

/* ===========================
 * Backend API contracts (NEW)
 * ===========================
 * These interfaces define the request/response shapes for the
 * ViennaRNA-backed stems analysis endpoint.
 */

export interface StemsRequest {
  /** DNA sequence (A/C/G/T preferred; case-insensitive) */
  sequence: string;
  /** Minimal final merged stem length to keep (bp) */
  min_stem_len: number;
  /** Merge adjacent stems when gap (unpaired) ≤ this value */
  merge_max_gap: number;
}

export interface StemsResponse {
  /** Length of the provided sequence */
  length: number;
  /** Merged stem regions suitable for direct visualization */
  regions: Region[];
}
