// File: frontend/src/app/features/feature-picker/feature-picker.component.ts
// Version: v2.1.0
//
// Compact “Detected features” panel.
// - Small chips with counts
// - Tiny status row for stems (loading / error)
// - Shorter scroll area for the ranges table
// - Same preferred order, includes 'stems'
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

  /** Preferred order; any extra kinds (future) are appended after */
  private readonly ORDER: FeatureKind[] = [
    'homopolymers',
    'gc-low',
    'gc-high',
    'entropy-low',
    'repeats',
    'long-repeats',
    'stems',
  ];

  kinds = computed<FeatureKind[]>(() => {
    const by = this.a.resultByKind();
    const keys = Object.keys(by) as FeatureKind[];
    const ordered = this.ORDER.filter(k => keys.includes(k));
    for (const k of keys) if (!this.ORDER.includes(k)) ordered.push(k);
    return ordered;
  });

  counts = computed<Record<FeatureKind, number>>(() => {
    const by = this.a.resultByKind();
    const out = {} as Record<FeatureKind, number>;
    for (const k of this.kinds()) out[k] = by[k]?.length ?? 0;
    return out;
  });

  ranges = computed<FeatureRegion[]>(() => this.a.selectedRegions());

  selectedIdx = computed<number>(() => {
    const sel = this.a.selectedRange();
    if (!sel) return -1;
    const list = this.a.selectedRegions();
    return list.findIndex(r => r.start === sel.start && r.end === sel.end && r.kind === sel.kind);
  });

  selected = computed<FeatureRegion | null>(() => this.a.selectedRange());
  kindColor = computed<string>(() => {
    const k = this.a.selectedKind();
    return k ? this.a.featureColor(k) : '#94a3b8';
  });

  selectKind(k: FeatureKind) {
    this.a.selectKind(k);
    this.a.selectRangeByIndex(null);
    if (k === 'stems') { void this.a.refreshStems(); }
  }
  clearKind() { this.a.selectKind(null); }
  pick(i: number) { this.a.selectRangeByIndex(i); }
  clearRange() { this.a.selectRange(null); }

  len(r: FeatureRegion) { return r.end - r.start; }
}
