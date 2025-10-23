// File: frontend/src/app/primers/models/primer-params.model.ts
// Version: v0.1.0
/**
 * Models for primer parameters and request/response DTOs.
 * Mirrors backend Pydantic schemas for tight typing.
 */

export interface PrimerDesignWeights {
  wTm: number;
  wGC: number;
  wDimer: number;
  w3p: number;
  wOff: number;
}

export interface PrimerDesignParameters {
  primerLengthMin: number;
  primerLengthMax: number;
  primerTmMin: number;
  primerTmMax: number;
  primerGCMin: number;
  primerGCMax: number;
  primerHomopolymerMax: number;
  primerThreePrimeEndLength: number;
  primerThreePrimePrimerMatchMax: number;
  primerTargetMatchMax: number;
  primerTargetMatchNumberMax: number;
  primerSecondaryStructureDeltaGMin: number;
  primerTmDifferenceMax: number;
  weights: PrimerDesignWeights;
  randomSeed?: number | null;
}

export interface TargetAnalysisParameters {
  targetGCMin: number;
  targetGCMax: number;
  targetSlidingWindowSize: number;
  targetHomopolymerMax: number;
  targetRepeatMax: number;
}

export interface PrimerDesignRequest {
  sequence: string;
  seqTargetStartPosition: number; // 1-based inclusive
  seqTargetEndPosition: number;   // 1-based inclusive
  parameters: PrimerDesignParameters;
}

export interface PrimerSeqInfo {
  sequence: string;
  tm: number;
  gc: number;
}

export interface PrimerDesignResponse {
  runId: string;
  forwardPrimer: PrimerSeqInfo;
  reversePrimer: PrimerSeqInfo;
  pairScore: number;
  warnings: string[];
}

export interface TargetAnalysisRequest {
  sequence: string;
  parameters: TargetAnalysisParameters;
}

export interface TargetAnalysisResponse {
  gcDistribution: number[];
  homopolymers: Array<Record<string, unknown>>;
  repeats: Array<Record<string, unknown>>;
  structureWarnings: string[];
}

export interface AlignmentsRequest {
  forward: string;
  reverse: string;
  target: string;
}

export interface AlignmentsResponse {
  forwardTarget: string;
  reverseTarget: string;
  selfDimer: string;
  crossDimer: string;
}
