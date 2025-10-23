// File: frontend/src/app/primers/models/primer-run.model.ts
// Version: v0.1.0
/**
 * Models for run records used in history and report pages.
 */

import { PrimerDesignResponse } from './primer-params.model';

export interface PrimerRunRecord {
  id: string;
  createdAt: string;
  targetStart: number;
  targetEnd: number;
  sequenceDigest: string;
  status: 'completed' | 'failed' | string;
  result?: PrimerDesignResponse | null;
}
