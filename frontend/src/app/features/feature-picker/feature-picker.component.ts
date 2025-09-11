// File: frontend/src/app/features/feature-picker/feature-picker.component.ts
// Version: v2.0.0
//
// Replaces legacy Results UI:
// - Chip selector with counts
// - Pinned selected row (sticky)
// - Compact 3-column table (start, end, length) with vertical scroll only
// - Includes dynamic kinds with preferred order, including 'stems'
import { Component, computed } from '@angular/core';
import { CommonModule } from '@angular/common';
import { AnalysisService, FeatureKind, FeatureRegion } from '../../core/services/analysis.service';

@Component({
  selector: 'app-feature-picker',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './feature-picker.component.html',
  styleUrls: ['./feature-picker.component.css'],
})
export class FeaturePickerComponent {
  constructor(public a: AnalysisService) {}

  /** Preferred order; any extra kinds (future) are appended at the end */
  private readonly ORDER: FeatureKind[] = [
    'homopolymers',
    'gc-low',
    'gc-high',
    'entropy-low',
    'repeats',
    'long-repeats',
    'stems',
  ];

  /** Kinds to render (dynamic from service map, ordered by ORDER then any extras) */
  kinds = computed<FeatureKind[]>(() => {
    const map = this.a.resultByKind();
    const keys = Object.keys(map) as FeatureKind[];
    const seen = new Set(keys);
    const ordered = this.ORDER.filter(k => seen.has(k));
    for (const k of keys) if (!this.ORDER.includes(k)) ordered.push(k);
    return ordered;
  });

  /** Per-kind counts */
  counts = computed<Record<FeatureKind, number>>(() => {
    const by = this.a.resultByKind();
    const out = {} as Record<FeatureKind, number>;
    for (const k of this.kinds()) out[k] = by[k]?.length ?? 0;
    return out;
  });

  /** Ranges for the selected kind */
  ranges = computed<FeatureRegion[]>(() => this.a.selectedRegions());

  /** Selected index in the current list */
  selectedIdx = computed<number>(() => {
    const sel = this.a.selectedRange();
    if (!sel) return -1;
    const list = this.a.selectedRegions();
    return list.findIndex(r => r.start === sel.start && r.end === sel.end && r.kind === sel.kind);
  });

  /** Convenience */
  selected = computed<FeatureRegion | null>(() => this.a.selectedRange());
  kindColor = computed<string>(() => {
    const k = this.a.selectedKind();
    return k ? this.a.featureColor(k) : '#94a3b8';
  });

  /** Actions */
  selectKind(k: FeatureKind) {
    this.a.selectKind(k);
    this.a.selectRangeByIndex(null);
    if (k === 'stems') { void this.a.refreshStems(); } // ensure up-to-date from server
  }
  clearKind() { this.a.selectKind(null); }
  pick(i: number) { this.a.selectRangeByIndex(i); }
  clearRange() { this.a.selectRange(null); }

  /** Display helper: compute 0-based length */
  len(r: FeatureRegion) { return r.end - r.start; }
}
