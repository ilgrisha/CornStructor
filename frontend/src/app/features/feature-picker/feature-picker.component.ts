/**
 * Path: frontend/src/app/features/feature-picker/feature-picker.component.ts
 * Version: v1.0.0
 *
 * Feature/region chooser:
 *  - Click a feature type → loads its detected regions
 *  - Click a region → highlights it in the sequence view
 */
import { Component, computed } from '@angular/core';
import { NgFor, NgIf, NgClass } from '@angular/common';
import { AnalysisService, FeatureKind } from '../../core/services/analysis.service';

type Row = { start: number; end: number; len: number };

@Component({
  selector: 'app-feature-picker',
  standalone: true,
  imports: [NgFor, NgIf, NgClass],
  templateUrl: './feature-picker.component.html',
  styleUrls: ['./feature-picker.component.css'],
})
export class FeaturePickerComponent {
  constructor(public a: AnalysisService) {}

  kinds: FeatureKind[] = ['invalid','homopolymers','gcLow','gcHigh','entropyLow','repeats'];

  counts = computed(() => {
    const f = this.a.features();
    return {
      invalid: f.invalid.length,
      homopolymers: f.homopolymers.length,
      gcLow: f.gcLow.length,
      gcHigh: f.gcHigh.length,
      entropyLow: f.entropyLow.length,
      repeats: f.repeats.length,
    } as Record<FeatureKind, number>;
  });

  rows = computed<Row[]>(() => {
    const k = this.a.selectedFeature();
    if (!k) return [];
    return this.a.getFeatureArray(k).map(r => ({ start: r.start, end: r.end, len: r.end - r.start }));
  });

  pickKind(k: FeatureKind) { this.a.selectFeature(k); }
  pickRow(i: number) { this.a.selectRegion(i); }

  isKindActive(k: FeatureKind) { return this.a.selectedFeature() === k; }
  isRowActive(i: number) { return this.a.selectedIndex() === i; }

  color(k: FeatureKind) { return this.a.featureColors[k]; }
}
