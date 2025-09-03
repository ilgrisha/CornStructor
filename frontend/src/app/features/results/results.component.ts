// ==AI-HEADER-START==
// Path: frontend/src/app/features/results/results.component.ts
// Version: 1
// Patch: 20250903-065943-a90b43
// Stamp: 2025-09-03T06:59:43Z
// ==AI-HEADER-END==
/* ============================================================================

Path: frontend/src/app/features/results/results.component.ts

Version: v2.4.0

============================================================================

Compact 3-column table (start, end, length) with no horizontal scroll

Table has its own vertical scroll; card itself does not scroll

Selected range is pinned at the top of the table (sticky), so it doesn't

move while the user scrolls the list

Removed scrollIntoView to avoid jumping on selection

==========================================================================*/
import { Component, computed } from '@angular/core';
import { CommonModule } from '@angular/common';
import { AnalysisService, FeatureKind, FeatureRegion } from '../../core/services/analysis.service';

@Component({
selector: 'app-results',
standalone: true,
imports: [CommonModule],
templateUrl: './results.component.html',
styleUrls: ['./results.component.css'],
})
export class ResultsComponent {
constructor(public a: AnalysisService) {}

kinds: FeatureKind[] = ['homopolymers', 'gc-low', 'gc-high', 'entropy-low', 'repeats', 'long-repeats'];

counts = computed(() => {
const r = this.a.resultByKind();
return {
'homopolymers': r['homopolymers'].length,
'gc-low': r['gc-low'].length,
'gc-high': r['gc-high'].length,
'entropy-low': r['entropy-low'].length,
'repeats': r['repeats'].length,
'long-repeats': r['long-repeats'].length,
} as Record<FeatureKind, number>;
});

ranges = computed<FeatureRegion[]>(() => this.a.selectedRegions());

selectedIdx = computed(() => {
const sel = this.a.selectedRange();
if (!sel) return -1;
const list = this.a.selectedRegions();
return list.findIndex(r => r.start === sel.start && r.end === sel.end);
});

selected = computed<FeatureRegion | null>(() => this.a.selectedRange());
kindColor = computed(() => {
const k = this.a.selectedKind();
return k ? this.a.featureColor(k) : '#94a3b8';
});

selectKind(k: FeatureKind) { this.a.selectKind(k); }
clearKind() { this.a.selectKind(null); }

pick(idx: number) { this.a.selectRangeByIndex(idx); }
clearRange() { this.a.selectRange(null); }

len(r: FeatureRegion) { return r.end - r.start; } // 0-based [start,end)
}
