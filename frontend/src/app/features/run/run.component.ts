/* ============================================================================
 * Path: frontend/src/app/features/run/run.component.ts
 * Version: v2.1.0
 * ==========================================================================*/
import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { ApiService } from '../../core/services/api.service';
import { AnalysisService } from '../../core/services/analysis.service';

@Component({
  selector: 'app-run',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './run.component.html',
  styleUrls: ['./run.component.css']
})
export class RunComponent {
  constructor(public api: ApiService, private a: AnalysisService) {}

  run() {
    const req = {
      sequence: this.a.sequence(),
      params: this.a.params(),
      toggles: { homopolymers: true, gc: true, entropy: true, repeats: true }
    };
    this.api.startDesign(req);
  }
}
