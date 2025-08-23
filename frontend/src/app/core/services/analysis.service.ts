/* ============================================================================
 * Path: frontend/src/app/core/services/analysis.service.ts
 * Version: v3.1.0
 * ============================================================================
 * Holds sequence + caret (signals), analysis params/toggles, detectors, and
 * selection API. Now includes alias param names `gcWindow`/`entropyWindow`
 * (kept in sync with `gcWin`/`entropyWin`). Adds `setSequence(...)`.
 * ==========================================================================*/

import { Injectable, WritableSignal, Signal, computed, signal } from '@angular/core';

export type FeatureKind =
  | 'invalid'
  | 'homopolymers'
  | 'gcLow'
  | 'gcHigh'
  | 'entropyLow'
  | 'repeats';

export interface Interval {
  start: number; // inclusive
  end: number;   // exclusive
  meta?: Record<string, unknown>;
}

export interface AnalysisParams {
  wrap: number;

  homoMin: number;

  // GC window aliases (templates can use either; both exist on the object)
  gcWin: number;
  gcWindow: number;
  gcMin: number;
  gcMax: number;

  // Entropy window aliases
  entropyWin: number;
  entropyWindow: number;
  entropyMin: number;

  repeatKmin: number;
  repeatKmax: number;
  repeatMinReps: number;
}

export interface Toggles {
  homopolymers: boolean;
  gc: boolean;
  entropy: boolean;
  repeats: boolean;
}

@Injectable({ providedIn: 'root' })
export class AnalysisService {
  // --- Sequence & caret ------------------------------------------------------
  sequence: WritableSignal<string> = signal('');
  caret: WritableSignal<number | null> = signal<number | null>(null);
  length = computed(() => this.sequence().length);

  setCaret(pos: number | null) { this.caret.set(pos); }

  /** Public setter used by sequence-input component */
  setSequence(seq: string) {
    // Keep uppercase consistently
    this.sequence.set((seq || '').toUpperCase());
  }

  // --- Parameters & toggles --------------------------------------------------
  params: WritableSignal<AnalysisParams> = signal<AnalysisParams>({
    wrap: 60,

    homoMin: 5,

    gcWin: 100,
    gcWindow: 100,
    gcMin: 40,
    gcMax: 60,

    entropyWin: 100,
    entropyWindow: 100,
    entropyMin: 1.6, // bits (0..~2 for DNA)

    repeatKmin: 2,
    repeatKmax: 6,
    repeatMinReps: 3,
  });

  toggles: WritableSignal<Toggles> = signal<Toggles>({
    homopolymers: true,
    gc: true,
    entropy: true,
    repeats: true,
  });

  updateParam<K extends keyof AnalysisParams>(key: K, value: AnalysisParams[K]) {
    this.params.update(p => {
      const next: AnalysisParams = { ...p, [key]: value } as AnalysisParams;

      // Keep window name aliases in sync
      if (key === 'gcWin')        next.gcWindow = Number(value);
      if (key === 'gcWindow')     next.gcWin    = Number(value);
      if (key === 'entropyWin')   next.entropyWindow = Number(value);
      if (key === 'entropyWindow')next.entropyWin    = Number(value);

      return next;
    });
  }

  updateToggle<K extends keyof Toggles>(key: K, value: Toggles[K]) {
    this.toggles.update(t => ({ ...t, [key]: value }));
  }

  // --- Feature colors --------------------------------------------------------
  readonly featureColors: Record<FeatureKind, string> = {
    invalid:     '#FF3B30',
    homopolymers:'#FFD60A',
    gcLow:       '#64D2FF',
    gcHigh:      '#FF9FF3',
    entropyLow:  '#FF9F0A',
    repeats:     '#34C759',
  };

  // --- Selection state -------------------------------------------------------
  selectedFeature: WritableSignal<FeatureKind | null> = signal<FeatureKind | null>(null);
  selectedIndex:  WritableSignal<number | null> = signal<number | null>(null);

  selectFeature(kind: FeatureKind | null) {
    this.selectedFeature.set(kind);
    this.selectedIndex.set(null);
  }
  selectRegion(idx: number | null) {
    if (this.selectedFeature() == null) { this.selectedIndex.set(null); return; }
    this.selectedIndex.set(idx);
  }
  currentRegion = computed<({ kind: FeatureKind; start: number; end: number } | null)>(() => {
    const k = this.selectedFeature();
    const i = this.selectedIndex();
    if (k == null || i == null) return null;
    const arr = this.getFeatureArray(k);
    const it = arr[i];
    return it ? { kind: k, start: it.start, end: it.end } : null;
  });

  // --- Analysis outputs ------------------------------------------------------
  invalid: Signal<Interval[]> = computed(() => {
    const s = this.sequence().toUpperCase();
    const bad: Interval[] = [];
    let i = 0;
    while (i < s.length) {
      if (!'ACGT'.includes(s[i])) {
        const start = i;
        while (i < s.length && !'ACGT'.includes(s[i])) i++;
        bad.push({ start, end: i });
      } else {
        i++;
      }
    }
    return bad;
  });

  homopolymers: Signal<Interval[]> = computed(() => {
    const s = this.sequence().toUpperCase();
    const min = Math.max(2, this.params().homoMin);
    const out: Interval[] = [];
    let i = 0;
    while (i < s.length) {
      const ch = s[i];
      if (!'ACGT'.includes(ch)) { i++; continue; }
      let j = i + 1;
      while (j < s.length && s[j] === ch) j++;
      if (j - i >= min) out.push({ start: i, end: j, meta: { base: ch, len: j - i } });
      i = j;
    }
    return out;
  });

  gcLow: Signal<Interval[]> = computed(() => this._gcOut('low'));
  gcHigh: Signal<Interval[]> = computed(() => this._gcOut('high'));

  entropyLow: Signal<Interval[]> = computed(() => {
    const { entropyWin, entropyMin } = this.params();
    const s = this.sequence().toUpperCase();
    const L = s.length;
    if (!L) return [];
    const win = Math.max(1, Math.min(entropyWin, L));

    let a=0,c=0,g=0,t=0,n=0;
    const push = (ch: string, d: number) => {
      switch (ch) {
        case 'A': a+=d; break; case 'C': c+=d; break;
        case 'G': g+=d; break; case 'T': t+=d; break;
        default: n+=d;
      }
    };

    for (let i=0;i<win;i++) push(s[i], +1);

    const entropy = () => {
      const tot = a+c+g+t;
      if (tot === 0) return 0;
      const probs = [a,c,g,t].filter(x=>x>0).map(x=>x/tot);
      let H = 0;
      for (const p of probs) H += -p * Math.log2(p);
      return H;
    };

    const windows: Interval[] = [];
    if (entropy() < entropyMin) windows.push({ start: 0, end: win });

    for (let i=win; i<L; i++) {
      push(s[i-win], -1);
      push(s[i], +1);
      if (entropy() < entropyMin) {
        windows.push({ start: i - win + 1, end: i + 1 });
      }
    }
    return this._mergeIntervals(windows);
  });

  repeats: Signal<Interval[]> = computed(() => {
    const { repeatKmin, repeatKmax, repeatMinReps } = this.params();
    const s = this.sequence().toUpperCase();
    const out: Interval[] = [];
    const L = s.length;
    if (!L) return out;

    for (let pos = 0; pos < L; pos++) {
      for (let k = repeatKmin; k <= repeatKmax; k++) {
        if (pos + k * repeatMinReps > L) break;
        const unit = s.slice(pos, pos + k);
        if (!/^[ACGT]+$/.test(unit)) continue;
        let reps = 1;
        while (pos + k * (reps + 1) <= L &&
               s.slice(pos + k*reps, pos + k*(reps+1)) === unit) {
          reps++;
        }
        if (reps >= repeatMinReps) {
          out.push({ start: pos, end: pos + k*reps, meta: { k, reps, unit } });
          pos = pos + k*reps - 1;
          break;
        }
      }
    }
    return this._mergeIntervals(out);
  });

  features = computed(() => ({
    invalid: this.invalid(),
    homopolymers: this.homopolymers(),
    gcLow: this.gcLow(),
    gcHigh: this.gcHigh(),
    entropyLow: this.entropyLow(),
    repeats: this.repeats(),
  }));

  getFeatureArray(kind: FeatureKind): Interval[] {
    switch (kind) {
      case 'invalid':      return this.invalid();
      case 'homopolymers': return this.homopolymers();
      case 'gcLow':        return this.gcLow();
      case 'gcHigh':       return this.gcHigh();
      case 'entropyLow':   return this.entropyLow();
      case 'repeats':      return this.repeats();
    }
  }

  // --- Private helpers -------------------------------------------------------
  private _gcOut(which: 'low'|'high'): Interval[] {
    const { gcWin, gcMin, gcMax } = this.params();
    const s = this.sequence().toUpperCase();
    const L = s.length;
    if (!L) return [];

    const win = Math.max(1, Math.min(gcWin, L));
    const isGC = (ch: string) => ch === 'G' || ch === 'C';

    let gc = 0, valid = 0;
    for (let i=0;i<win;i++) { if ('ACGT'.includes(s[i])) { valid++; if (isGC(s[i])) gc++; } }
    const windows: Interval[] = [];
    const pushIf = (start: number) => {
      if (valid === 0) return;
      const pct = (gc/valid) * 100;
      if (which === 'low'  && pct < gcMin) windows.push({ start, end: start + win, meta: { pct } });
      if (which === 'high' && pct > gcMax) windows.push({ start, end: start + win, meta: { pct } });
    };
    pushIf(0);

    for (let i=win; i<L; i++) {
      const out = s[i-win], inn = s[i];
      if ('ACGT'.includes(out)) { valid--; if (isGC(out)) gc--; }
      if ('ACGT'.includes(inn)) { valid++; if (isGC(inn)) gc++; }
      pushIf(i - win + 1);
    }
    return this._mergeIntervals(windows);
  }

  private _mergeIntervals(arr: Interval[]): Interval[] {
    if (arr.length <= 1) return arr.slice();
    const sorted = arr.slice().sort((a,b) => a.start - b.start || a.end - b.end);
    const out: Interval[] = [];
    let cur = { ...sorted[0] };
    for (let i = 1; i < sorted.length; i++) {
      const it = sorted[i];
      if (it.start <= cur.end) {
        cur.end = Math.max(cur.end, it.end);
      } else {
        out.push(cur);
        cur = { ...it };
      }
    }
    out.push(cur);
    return out;
  }
}
