/* ============================================================================
 * Path: frontend/src/app/features/params/params.component.ts
 * Version: v2.3.1
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
  constructor(private a: AnalysisService) {}
  p = computed(() => this.a.params());
  updateNum(key: keyof AnalysisParams, ev: Event) {
    const v = Number((ev.target as HTMLInputElement).value);
    this.a.updateParam(key, Number.isFinite(v) ? v : this.p()[key]);
  }
}
