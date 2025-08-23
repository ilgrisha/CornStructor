export type FeatureKey = 'homopolymers' | 'gc' | 'entropy' | 'repeats';

export interface Region {
  start: number;   // 0-based inclusive
  end: number;     // 0-based exclusive
  label?: string;
  extra?: Record<string, unknown>;
}

export interface AnalysisParams {
  // homopolymers
  homoMin: number;

  // GC%
  gcWindow: number;
  gcMin: number;   // percent 0..100
  gcMax: number;

  // low complexity via Shannon entropy
  entropyWindow: number;
  entropyMin: number; // bits 0..2 (on A,C,G,T only)

  // simple tandem repeats
  repeatKmin: number;
  repeatKmax: number;
  repeatMinReps: number;
}

export interface FeatureToggles {
  homopolymers: boolean;
  gc: boolean;
  entropy: boolean;
  repeats: boolean;
}

export interface AnalysisResult {
  homopolymers: Region[];
  gcOutOfRange: Region[];
  lowComplexity: Region[];
  repeats: Region[];
}
