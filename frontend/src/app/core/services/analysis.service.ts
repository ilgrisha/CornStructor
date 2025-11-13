/* File: frontend/src/app/core/services/analysis.service.ts
/* Version: v0.5.0
/*==========================================================================
Stems server-call only on sequence changes; params re-filter locally.

- Calls ViennaRNA backend *only* when the sequence actually changes.
- Caches the "base" stems runs from server using (min_stem_len=1, merge_max_gap=0).
  These are raw consecutive paired runs derived from MFE dot-bracket.
- Applies UI parameters (stemMinLen, stemMergeMaxGap) *locally* to the cached runs,
  so changing parameters is instant and does NOT call the server.
- Keeps debounced auto-refresh on sequence edits/uploads (default 300 ms).
- Exposes stemsLoading() and stemsError() so UI can show status banners.

Other analyses (homopolymers, GC, entropy, repeats, long-repeats) unchanged.
==========================================================================*/
import { Injectable, computed, signal, WritableSignal, effect } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { firstValueFrom } from 'rxjs';

export type FeatureKind =
  | 'homopolymers'
  | 'gc-low'
  | 'gc-high'
  | 'entropy-low'
  | 'repeats'
  | 'long-repeats'
  | 'stems';

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

  longRepMinLen: number;
  longRepMinPct: number;

  // Secondary structure (client-side filtering of cached raw stems)
  stemMinLen: number;
  stemMergeMaxGap: number;
}

const API_BASE = '/api';

@Injectable({ providedIn: 'root' })
export class AnalysisService {
  constructor(private http: HttpClient) {}

  // ---- sequence state ----
  sequence: WritableSignal<string> = signal('');

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

  // ---- parameters ----
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
    longRepMinLen: 30,
    longRepMinPct: 85,
    stemMinLen: 4,
    stemMergeMaxGap: 2,
  });

  updateParam<K extends keyof AnalysisParams>(key: K, value: AnalysisParams[K]) {
    this.params.update(p => ({ ...p, [key]: value }));
  }

  featureColor(k: FeatureKind): string {
    switch (k) {
      case 'homopolymers': return '#f59e0b';
      case 'gc-low': return '#3b82f6';
      case 'gc-high': return '#10b981';
      case 'entropy-low': return '#ef4444';
      case 'repeats': return '#a855f7';
      case 'long-repeats': return '#f97316';
      case 'stems': return '#4b7bec';
      default: return '#64748b';
    }
  }

  // ---- local analyses (unchanged) ----
  private homopolymers = computed<FeatureRegion[]>(() => {
    const s = this.sequence().toUpperCase();
    const min = Math.max(2, Math.floor(this.params().homoMin));
    const out: FeatureRegion[] = [];
    let i = 0;
    while (i < s.length) {
      const ch = s[i];
      if (!(ch === 'A' || ch === 'C' || ch === 'G' || ch === 'T')) { i++; continue; }
      let j = i + 1;
      while (j < s.length && s[j] === ch) j++;
      if (j - i >= min) out.push({ kind: 'homopolymers', start: i, end: j, meta: { base: ch, len: j - i } });
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
      if (last && last.kind === kind && last.end >= start) last.end = Math.max(last.end, end);
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
    const win = Math.max(2, Math.floor(p.entropyWin));
    const thr = p.entropyMin;
    const out: FeatureRegion[] = [];
    if (!s.length || win > s.length) return out;

    const counts = { A: 0, C: 0, G: 0, T: 0 } as Record<string, number>;
    for (let i = 0; i < win; i++) counts[s[i]] = (counts[s[i]] ?? 0) + 1;

    const H = () => {
      const n = win;
      let h = 0;
      for (const b of ['A', 'C', 'G', 'T'] as const) {
        const p = (counts[b] || 0) / n;
        if (p > 0) h -= p * Math.log2(p);
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
      if (H() < thr) push(i, i + win);
      if (i + win < s.length) {
        counts[s[i]] = (counts[s[i]] || 0) - 1;
        counts[s[i + win]] = (counts[s[i + win]] || 0) + 1;
      }
    }
    return this.merge(out);
  });

  private repeats = computed<FeatureRegion[]>(() => {
    const s = this.sequence().toUpperCase();
    const p = this.params();
    const kmin = Math.max(2, Math.floor(p.repeatKmin));
    const kmax = Math.max(kmin, Math.floor(p.repeatKmax));
    const minReps = Math.max(2, Math.floor(p.repeatMinReps));
    const out: FeatureRegion[] = [];

    for (let i = 0; i < s.length; i++) {
      for (let k = kmin; k <= kmax && i + k <= s.length; k++) {
        const motif = s.slice(i, i + k);
        if (motif.split('').every(ch => ch === motif[0])) continue;
        let reps = 1;
        let j = i + k;
        while (j + k <= s.length && s.slice(j, j + k) === motif) {
          reps++; j += k;
        }
        if (reps >= minReps) {
          out.push({ kind: 'repeats', start: i, end: i + reps * k, meta: { motif, reps } });
          i = j - 1;
          break;
        }
      }
    }
    return this.merge(out);
  });

  private longRepeats = computed<FeatureRegion[]>(() => {
    const s = this.sequence().toUpperCase();
    const p = this.params();
    const Lmin = Math.max(10, Math.floor(p.longRepMinLen || 0));
    const minId = Math.min(100, Math.max(0, p.longRepMinPct || 0)) / 100;
    if (!s.length || s.length < Lmin) return [];

    const seed = Math.max(6, Math.min(12, Math.floor(Lmin / 3)));
    const index = new Map<string, number[]>();

    for (let i = 0; i + seed <= s.length; i++) {
      const kmer = s.slice(i, i + seed);
      if (kmer.split('').every(ch => ch === kmer[0])) continue;
      const arr = index.get(kmer);
      if (arr) { arr.push(i); if (arr.length > 200) index.set(kmer, arr.slice(0, 200)); }
      else index.set(kmer, [i]);
    }

    const out: FeatureRegion[] = [];
    const seen = new Set<string>();
    const maxCandidatesPerI = 25;
    const allowed = (len: number) => Math.floor(len * (1 - minId) + 1e-9);

    function extend(i: number, j: number) {
      let a = i, b = j;
      let len = 0, mism = 0;

      while (a + len < s.length && b + len < s.length) {
        if (s[a + len] !== s[b + len]) mism++;
        len++;
        if (mism > allowed(len)) {
          mism -= (s[a + len - 1] !== s[b + len - 1]) ? 1 : 0;
          len--;
          break;
        }
      }
      if (len < Lmin) return null;

      let la = a - 1, lb = b - 1;
      while (la >= 0 && lb >= 0) {
        const addMismatch = s[la] === s[lb] ? 0 : 1;
        const newLen = len + 1;
        if (mism + addMismatch <= allowed(newLen)) {
          a = la; b = lb; len = newLen; mism += addMismatch;
          la--; lb--;
        } else break;
      }
      const id = (len - mism) / len;
      return { a, b, len, id };
    }

    for (const pos of index.values()) {
      if (pos.length < 2) continue;
      for (let ii = 0; ii < pos.length; ii++) {
        const i = pos[ii];
        let tried = 0;
        for (let jj = ii + 1; jj < pos.length; jj++) {
          const j = pos[jj];
          if (j <= i) continue;
          const key = i + ':' + j;
          if (seen.has(key)) continue;
          const ex = extend(i, j);
          seen.add(key);
          if (ex) {
            const { a, b, len, id } = ex;
            if (len >= Lmin && id >= minId) {
              const pairId = a + '-' + b + '-' + len;
              out.push({ kind: 'long-repeats', start: a, end: a + len, meta: { mateStart: b, mateEnd: b + len, identity: +(id * 100).toFixed(1), pairId } });
              out.push({ kind: 'long-repeats', start: b, end: b + len, meta: { mateStart: a, mateEnd: a + len, identity: +(id * 100).toFixed(1), pairId } });
            }
          }
          tried++;
          if (tried >= maxCandidatesPerI) break;
        }
      }
    }

    return this.merge(out);
  });

  // ---- stems (server only on sequence change) ----
  /** Raw consecutive paired runs from the server (computed with min=1, gap=0), cached per-sequence. */
  private stemsBase = signal<FeatureRegion[]>([]);
  /** Loading & error for UI feedback */
  private _stemsLoading = signal<boolean>(false);
  private _stemsError = signal<string | null>(null);

  /** Expose status to components */
  stemsLoading = () => this._stemsLoading();
  stemsError = () => this._stemsError();

  /** Apply current UI params to cached raw stems locally (merge + min length). */
  private stemsFiltered = computed<FeatureRegion[]>(() => {
    const base = this.stemsBase();
    const { stemMinLen, stemMergeMaxGap } = this.params();

    // Merge adjacent if gap <= stemMergeMaxGap (kind is 'stems' for all).
    const merged = this.mergeWithGap(base, stemMergeMaxGap);
    // Filter by minimal length
    const kept = merged.filter(r => (r.end - r.start) >= Math.max(1, Math.floor(stemMinLen)));
    return kept;
  });

  /** Public map consumed by the UI */
  resultByKind = computed<Record<FeatureKind, FeatureRegion[]>>(() => ({
    homopolymers: this.homopolymers(),
    'gc-low': this.gcBands().lows,
    'gc-high': this.gcBands().highs,
    'entropy-low': this.entropyLow(),
    repeats: this.repeats(),
    'long-repeats': this.longRepeats(),
    stems: this.stemsFiltered(),
  }));

  // ---- selection state ----
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

  /** Manually force a stems recomputation from the server (uses raw settings). */
  async refreshStems(): Promise<void> {
    const seq = (this.sequence() || '').toUpperCase().replace(/\s+/g, '');
    if (!seq) { this.stemsBase.set([]); return; }

    // Call backend with minimal filtering so we cache raw paired runs.
    const body = { sequence: seq, min_stem_len: 1, merge_max_gap: 0 };

    try {
      this._stemsLoading.set(true);
      this._stemsError.set(null);
      const resp = await firstValueFrom(
        this.http.post<{ length: number; regions: { start: number; end: number }[] }>(
          `${API_BASE}/analysis/stems`,
          body
        )
      );
      const regions: FeatureRegion[] = (resp?.regions || []).map(r => ({
        kind: 'stems', start: r.start, end: r.end,
      }));
      // No further merging/filtering here; stemsFiltered() applies UI params on the fly.
      this.stemsBase.set(this.sortNormalize(regions));
    } catch (err: any) {
      console.error('Stems analysis failed:', err);
      this._stemsError.set(err?.message ?? 'Unknown error');
      this.stemsBase.set([]);
    } finally {
      this._stemsLoading.set(false);
    }
  }

  // ---- debounced auto-refresh when *sequence* changes ONLY ----
  private _stemsDebounce: any = null;
  private _stemsEffect = effect(
    () => {
      const _seq = this.sequence(); // dependency: sequence only
      void _seq;

      if (this._stemsDebounce) clearTimeout(this._stemsDebounce);
      this._stemsDebounce = setTimeout(() => {
        void this.refreshStems();
      }, 300); // debounce 300ms after last change
    },
    { allowSignalWrites: true }
  );

  // ---- helpers ----
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

  /** Merge adjacent intervals if the gap (start - prev.end) ≤ maxGap. Assumes same kind. */
  private mergeWithGap(list: FeatureRegion[], maxGap: number): FeatureRegion[] {
    if (list.length <= 1) return list.slice();
    const arr = this.sortNormalize(list);
    const out: FeatureRegion[] = [];
    let cur = { ...arr[0] };
    for (let i = 1; i < arr.length; i++) {
      const r = arr[i];
      const gap = r.start - cur.end;
      if (gap <= Math.max(0, Math.floor(maxGap))) {
        cur.end = Math.max(cur.end, r.end);
      } else {
        out.push(cur);
        cur = { ...r };
      }
    }
    out.push(cur);
    return out;
  }

  private sortNormalize(list: FeatureRegion[]): FeatureRegion[] {
    const arr = list.slice().sort((a, b) => a.start - b.start || a.end - b.end);
    // Also clamp negatives and ensure start <= end
    return arr.map(r => ({
      kind: r.kind,
      start: Math.max(0, Math.min(r.start, r.end)),
      end: Math.max(0, Math.max(r.start, r.end)),
      meta: r.meta
    }));
  }
}
