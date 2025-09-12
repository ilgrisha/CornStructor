/* ============================================================================
 * Path: frontend/src/app/features/params/params.component.ts
 * Version: v2.4.0
 * ============================================================================
 * Fix: make AnalysisService instance PUBLIC so the template can call
 * `a.refreshStems()` without TS2341 ("property 'a' is private").
 * Also adds a tiny helper `refreshStemsNow()` if you prefer using a method
 * from the template instead of accessing `a` directly.
 * ==========================================================================*/
import { Component, computed } from '@angular/core';
import { CommonModule } from '@angular/common';
import { AnalysisService, AnalysisParams } from '../../core/services/analysis.service';

@Component({
  selector: 'app-params',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './params.component.html',
  styleUrls: ['./params.component.css'],
})
export class ParamsComponent {
  // NOTE: must be public for template access
  constructor(public a: AnalysisService) {}

  p = computed(() => this.a.params());

  updateNum(key: keyof AnalysisParams, ev: Event) {
    const v = Number((ev.target as HTMLInputElement).value);
    this.a.updateParam(key, Number.isFinite(v) ? v : this.p()[key]);
  }

  /** Optional helper if you'd rather call this from the template */
  refreshStemsNow() {
    void this.a.refreshStems();
  }
}
