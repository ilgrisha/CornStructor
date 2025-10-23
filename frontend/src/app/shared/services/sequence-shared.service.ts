// File: frontend/src/app/shared/services/sequence-shared.service.ts
// Version: v0.1.0
/**
 * SequenceSharedService
 * ---------------------
 * A simple cross-feature store for the currently selected sequence and target window.
 * This allows the Construction Tree tab and the Primers tab to share the Sequence box.
 */

import { Injectable } from '@angular/core';
import { BehaviorSubject } from 'rxjs';

export interface SharedSequenceState {
  sequence: string;          // raw sequence (no FASTA headers)
  start: number | null;      // 1-based inclusive
  end: number | null;        // 1-based inclusive
  name?: string | null;      // optional label (e.g., record id or FASTA header)
}

@Injectable({ providedIn: 'root' })
export class SequenceSharedService {
  private state$ = new BehaviorSubject<SharedSequenceState>({
    sequence: '',
    start: null,
    end: null,
    name: null
  });

  readonly value$ = this.state$.asObservable();

  set(state: Partial<SharedSequenceState>) {
    this.state$.next({ ...this.state$.value, ...state });
  }

  get snapshot(): SharedSequenceState {
    return this.state$.value;
  }

  clear() {
    this.state$.next({ sequence: '', start: null, end: null, name: null });
  }
}
