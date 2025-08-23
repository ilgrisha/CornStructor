/**
 * Path: frontend/src/app/features/results/results.component.ts
 * Version: v0.8.0
 *
 * Results card: places FeaturePicker (left) and SequenceView (right),
 * imports DecimalPipe so {{ len() | number }} works.
 */
import { Component, computed } from '@angular/core';
import { CardComponent } from '../../shared/ui/card/card.component';
import { NgIf, DecimalPipe } from '@angular/common';
import { AnalysisService } from '../../core/services/analysis.service';
import { SequenceViewComponent } from '../sequence-view/sequence-view.component';
import { FeaturePickerComponent } from '../feature-picker/feature-picker.component';

@Component({
  selector: 'app-results',
  standalone: true,
  imports: [CardComponent, NgIf, DecimalPipe, SequenceViewComponent, FeaturePickerComponent],
  templateUrl: './results.component.html',
  styleUrls: ['./results.component.css'],
})
export class ResultsComponent {
  constructor(public a: AnalysisService) {}
  len = this.a.length;
  stats = computed(() => {
    const f = this.a.features();
    return {
      invalid: f.invalid.length,
      homos: f.homopolymers.length,
      gcLow: f.gcLow.length,
      gcHigh: f.gcHigh.length,
      entropyLow: f.entropyLow.length,
      repeats: f.repeats.length,
    };
  });
}
