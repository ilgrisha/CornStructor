// ==AI-HEADER-START==
// Path: frontend/src/app/core/services/analysis.service.ts
// Version: 1
// Patch: 20250903-065943-a90b43
// Stamp: 2025-09-03T06:59:43Z
// ==AI-HEADER-END==
/* ============================================================================

Path: frontend/src/app/core/services/analysis.service.ts

Version: v3.4.0

============================================================================

Always compute all analyses

Separate GC LOW/HIGH; repeats exclude homopolymers

Overlap merging; selection signals (kind + single range)

Dynamic windows: gcWin, entropyWin

==========================================================================*/
import { Injectable, computed, signal, WritableSignal } from '@angular/core';

export type FeatureKind = 'homopolymers' | 'gc-low' | 'gc-high' | 'entropy-low' | 'repeats' | 'long-repeats';

export interface FeatureRegion {
kind: FeatureKind;
start: number; // 0-based inclusive
end: number; // 0-based exclusive
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

// Long repeats (approximate direct repeats)
longRepMinLen: number; // minimum aligned length (bp)
longRepMinPct: number; // identity percent 0..100
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

// params
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
});

updateParam<K extends keyof AnalysisParams>(key: K, value: AnalysisParams[K]) {
this.params.update(p => ({ ...p, [key]: value }));
}

// Color legend (UI)
featureColor(k: FeatureKind): string {
switch (k) {
case 'homopolymers': return '#f59e0b'; // amber
case 'gc-low': return '#3b82f6'; // blue
case 'gc-high': return '#10b981'; // emerald
case 'entropy-low': return '#ef4444'; // red
case 'repeats': return '#a855f7'; // purple
case 'long-repeats': return '#f97316'; // orange
default: return '#64748b';
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

// rolling base counts
const counts = { A: 0, C: 0, G: 0, T: 0 } as Record<string, number>;
for (let i = 0; i < win; i++) counts[s[i]] = (counts[s[i]] ?? 0) + 1;

const H = () => {
  const n = win;
  let h = 0;
  for (const b of ['A', 'C', 'G', 'T'] as const) {
    const p = counts[b] / n;
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
    counts[s[i]]--;
    counts[s[i + win]] = (counts[s[i + win]] ?? 0) + 1;
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
    if (motif.split('').every(ch => ch === motif[0])) continue; // skip homopolymer motifs
    let reps = 1;
    let j = i + k;
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

private longRepeats = computed<FeatureRegion[]>(() => {
// Approximate direct repeats (no gaps), report both occurrences as separate regions.
const s = this.sequence().toUpperCase();
const p = this.params();
const Lmin = Math.max(10, Math.floor(p.longRepMinLen || 0));
const minId = Math.min(100, Math.max(0, p.longRepMinPct || 0)) / 100;
if (!s.length || s.length < Lmin) return [];

const seed = Math.max(6, Math.min(12, Math.floor(Lmin / 3))); // heuristic seed length
const index = new Map<string, number[]>();

// Build seed index
for (let i = 0; i + seed <= s.length; i++) {
  const kmer = s.slice(i, i + seed);
  // skip pure homopolymer seeds to limit explosion
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

  // extend right to satisfy Lmin and as far as threshold allows
  while (a + len < s.length && b + len < s.length) {
    if (s[a + len] !== s[b + len]) mism++;
    len++;
    if (mism > allowed(len)) { // too many mismatches → back off last char
      mism -= (s[a + len - 1] !== s[b + len - 1]) ? 1 : 0;
      len--;
      break;
    }
  }
  // ensure at least Lmin
  if (len < Lmin) return null;

  // extend left while threshold holds
  let la = a - 1, lb = b - 1;
  while (la >= 0 && lb >= 0) {
    const addMismatch = s[la] === s[lb] ? 0 : 1;
    const newLen = len + 1;
    if (mism + addMismatch <= allowed(newLen)) {
      a = la; b = lb; len = newLen; mism += addMismatch;
      la--; lb--;
    } else break;
  }

  // final identity
  const id = (len - mism) / len;
  return { a, b, len, id };
}

// Iterate seeds
for (const [kmer, pos] of index) {
  if (pos.length < 2) continue;
  // For each left anchor, try a few right anchors
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
// Combined result by kind
resultByKind = computed<Record<FeatureKind, FeatureRegion[]>>(() => ({
homopolymers: this.homopolymers(),
'gc-low': this.gcBands().lows,
'gc-high': this.gcBands().highs,
'entropy-low': this.entropyLow(),
repeats: this.repeats(),
'long-repeats': this.longRepeats(),
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
