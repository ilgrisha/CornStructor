/**
 * Path: frontend/src/app/features/run/run.component.ts
 * Version: v0.5.0
 *
 * Bottom card: kicks off Construction Tree design and streams logs via SSE.
 * Uses ApiService.startDesign(request) and shows a link to the result page.
 */

import { Component } from '@angular/core';
import { CardComponent } from '../../shared/ui/card/card.component';
import { NgIf, NgFor } from '@angular/common';
import { ApiService } from '../../core/services/api.service';
import { AnalysisService } from '../../core/services/analysis.service';

@Component({
  selector: 'app-run',
  standalone: true,
  imports: [CardComponent, NgIf, NgFor],
  templateUrl: './run.component.html',
  styleUrls: ['./run.component.css'],
})
export class RunComponent {
  busy = false;

  constructor(public api: ApiService, public a: AnalysisService) {}

  start() {
    if (this.busy) return;
    this.busy = true;

    const req = {
      sequence: this.a.sequence(),
      params: this.a.params(),
      toggles: this.a.toggles(),
      // name: 'optional-tag'
    };

    this.api.startDesign(req)
      .catch(() => {
        const copy = [...this.api.logLines()];
        copy.push('Error starting design');
        this.api.logLines.set(copy);
      })
      .finally(() => (this.busy = false));
  }
}
