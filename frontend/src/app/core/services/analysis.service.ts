/* ============================================================================
 * Path: frontend/src/app/core/services/analysis.service.ts
 * Version: v3.3.0
 * ============================================================================
 * - Always compute all analyses
 * - Separate GC LOW/HIGH; repeats exclude homopolymers
 * - Overlap merging; selection signals (kind + single range)
 * - Dynamic windows: gcWin, entropyWin
 * ==========================================================================*/
import { Injectable, computed, signal, WritableSignal } from '@angular/core';

export type FeatureKind = 'homopolymers' | 'gc-low' | 'gc-high' | 'entropy-low' | 'repeats';

export interface FeatureRegion {
  kind: FeatureKind;
  start: number; // 0-based inclusive
  end: number;   // 0-based exclusive
  meta?: any;
}

export interface AnalysisParams {
  homoMin: number;

  gcWin: number;
  gcMin: number;
  gcMax: number;

  entropyWin: number;
  entropyMin: number;

  repeatKmin: number;
  repeatKmax: number;
  repeatMinReps: number;
}

@Injectable({ providedIn: 'root' })
export class AnalysisService {
  // raw sequence (uppercase enforced by views)
  sequence: WritableSignal<string> = signal('');

  // non-ACGT flags
  invalid = computed(() => {
    const s = this.sequence().toUpperCase();
    const out: boolean[] = new Array(s.length);
    for (let i = 0; i < s.length; i++) {
      const ch = s[i];
      out[i] = !(ch === 'A' || ch === 'C' || ch === 'G' || ch === 'T');
    }
    return out;
  });

  length = computed(() => this.sequence().length);

  params: WritableSignal<AnalysisParams> = signal({
    homoMin: 6,
    gcWin: 100,
    gcMin: 40,
    gcMax: 60,
    entropyWin: 50,
    entropyMin: 1.8,
    repeatKmin: 2,
    repeatKmax: 7,
    repeatMinReps: 3,
  });

  updateParam<K extends keyof AnalysisParams>(key: K, value: AnalysisParams[K]) {
    this.params.update(p => ({ ...p, [key]: value }));
  }

  // Color legend (UI)
  featureColor(k: FeatureKind): string {
    switch (k) {
      case 'homopolymers': return '#f59e0b';  // amber
      case 'gc-low':       return '#3b82f6';  // blue
      case 'gc-high':      return '#10b981';  // emerald
      case 'entropy-low':  return '#ef4444';  // red
      case 'repeats':      return '#a855f7';  // purple
      default:             return '#64748b';
    }
  }

  // --- Analyses (computed) ---
  private homopolymers = computed<FeatureRegion[]>(() => {
    const s = this.sequence().toUpperCase();
    const min = Math.max(2, Math.floor(this.params().homoMin));
    const out: FeatureRegion[] = [];
    let i = 0;
    while (i < s.length) {
      const ch = s[i];
      let j = i + 1;
      while (j < s.length && s[j] === ch) j++;
      const len = j - i;
      if ((ch === 'A' || ch === 'C' || ch === 'G' || ch === 'T') && len >= min) {
        out.push({ kind: 'homopolymers', start: i, end: j, meta: { nt: ch, len } });
      }
      i = j;
    }
    return this.merge(out);
  });

  private gcBands = computed<{ lows: FeatureRegion[]; highs: FeatureRegion[] }>(() => {
    const s = this.sequence().toUpperCase();
    const p = this.params();
    const win = Math.max(1, Math.floor(p.gcWin));
    const lows: FeatureRegion[] = [];
    const highs: FeatureRegion[] = [];
    if (!s.length || win > s.length) return { lows, highs };

    let gc = 0;
    for (let i = 0; i < win; i++) if (s[i] === 'G' || s[i] === 'C') gc++;
    const push = (arr: FeatureRegion[], kind: FeatureKind, start: number, end: number) => {
      if (end <= start) return;
      const last = arr[arr.length - 1];
      if (last && last.end >= start) last.end = Math.max(last.end, end);
      else arr.push({ kind, start, end });
    };

    for (let i = 0; i + win <= s.length; i++) {
      const pct = (gc / win) * 100;
      if (pct < p.gcMin) push(lows, 'gc-low', i, i + win);
      else if (pct > p.gcMax) push(highs, 'gc-high', i, i + win);

      if (i + win < s.length) {
        const out = s[i];
        const inn = s[i + win];
        if (out === 'G' || out === 'C') gc--;
        if (inn === 'G' || inn === 'C') gc++;
      }
    }
    return { lows: this.merge(lows), highs: this.merge(highs) };
  });

  private entropyLow = computed<FeatureRegion[]>(() => {
    const s = this.sequence().toUpperCase();
    const p = this.params();
    const win = Math.max(1, Math.floor(p.entropyWin));
    const out: FeatureRegion[] = [];
    if (!s.length || win > s.length) return out;

    const freq: Record<string, number> = { A: 0, C: 0, G: 0, T: 0 };
    const add = (ch: string, d: number) => { if (freq[ch] !== undefined) freq[ch] += d; };
    for (let i = 0; i < win; i++) add(s[i], 1);

    const H = () => {
      const tot = freq.A + freq.C + freq.G + freq.T;
      if (!tot) return 0;
      let h = 0;
      for (const k of ['A', 'C', 'G', 'T']) {
        const p = (freq as any)[k] / tot;
        if (p > 0) h += -p * Math.log2(p);
      }
      return h;
    };
    const push = (start: number, end: number) => {
      if (end <= start) return;
      const last = out[out.length - 1];
      if (last && last.end >= start) last.end = Math.max(last.end, end);
      else out.push({ kind: 'entropy-low', start, end });
    };

    for (let i = 0; i + win <= s.length; i++) {
      if (H() < p.entropyMin) push(i, i + win);
      if (i + win < s.length) {
        add(s[i], -1);
        add(s[i + win], +1);
      }
    }
    return this.merge(out);
  });

  private repeats = computed<FeatureRegion[]>(() => {
    // Tandem repeats, exclude homopolymers (k must be >=2)
    const s = this.sequence().toUpperCase();
    const p = this.params();
    const kmin = Math.max(2, Math.floor(p.repeatKmin));
    const kmax = Math.max(kmin, Math.floor(p.repeatKmax));
    const minReps = Math.max(2, Math.floor(p.repeatMinReps));
    const out: FeatureRegion[] = [];

    for (let i = 0; i < s.length; i++) {
      for (let k = kmin; k <= kmax && i + k <= s.length; k++) {
        const motif = s.slice(i, i + k);
        if (!/^[ACGT]+$/.test(motif)) continue;
        if (new Set(motif).size === 1) continue; // exclude homopolymers
        let j = i + k, reps = 1;
        while (j + k <= s.length && s.slice(j, j + k) === motif) {
          reps++; j += k;
        }
        if (reps >= minReps) {
          out.push({ kind: 'repeats', start: i, end: i + reps * k, meta: { motif, reps } });
          i = j - 1; // skip forward
          break;
        }
      }
    }
    return this.merge(out);
  });

  // Combined result by kind
  resultByKind = computed<Record<FeatureKind, FeatureRegion[]>>(() => ({
    homopolymers: this.homopolymers(),
    'gc-low': this.gcBands().lows,
    'gc-high': this.gcBands().highs,
    'entropy-low': this.entropyLow(),
    repeats: this.repeats(),
  }));

  // --- Selection state used across UI (sequence box + results table) ---
  selectedKind: WritableSignal<FeatureKind | null> = signal(null);
  selectedRange: WritableSignal<FeatureRegion | null> = signal(null);

  selectedRegions = computed<FeatureRegion[]>(() => {
    const k = this.selectedKind();
    return k ? this.resultByKind()[k] ?? [] : [];
  });

  selectKind(k: FeatureKind | null) {
    this.selectedKind.set(k);
    this.selectedRange.set(null);
  }
  selectRangeByIndex(idx: number | null) {
    if (idx === null) { this.selectedRange.set(null); return; }
    const list = this.selectedRegions();
    if (idx >= 0 && idx < list.length) this.selectedRange.set(list[idx]);
  }
  selectRange(r: FeatureRegion | null) {
    this.selectedRange.set(r);
  }

  // helpers
  private merge(list: FeatureRegion[]): FeatureRegion[] {
    if (list.length <= 1) return list.slice();
    const arr = list.slice().sort((a, b) => a.start - b.start);
    const out: FeatureRegion[] = [];
    let cur = { ...arr[0] };
    for (let i = 1; i < arr.length; i++) {
      const r = arr[i];
      if (r.start <= cur.end && r.kind === cur.kind) {
        cur.end = Math.max(cur.end, r.end);
      } else {
        out.push(cur);
        cur = { ...r };
      }
    }
    out.push(cur);
    return out;
  }
}
