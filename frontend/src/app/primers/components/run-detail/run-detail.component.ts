// File: frontend/src/app/primers/components/run-detail/run-detail.component.ts
// Version: v0.1.0
/**
 * RunDetailComponent (Angular Report)
 * -----------------------------------
 * Displays a single design run:
 *  - Primer sequences + Tm/GC + score + warnings
 *  - Alignments (forward↔target, reverse↔targetRC, self/cross-dimer)
 *  - Target analysis (GC windows and flags)
 * The "Report" is an Angular view, not a static HTML file.
 */

import { Component, OnInit } from '@angular/core';
import { ActivatedRoute } from '@angular/router';
import { PrimersService } from '../../services/primers.service';
import { PrimerRunRecord } from '../../models/primer-run.model';
import { AlignmentsResponse, TargetAnalysisRequest, TargetAnalysisResponse } from '../../models/primer-params.model';
import { SequenceSharedService } from '../../../shared/services/sequence-shared.service';

@Component({
  selector: 'app-run-detail',
  templateUrl: './run-detail.component.html',
  styleUrls: ['./run-detail.component.scss'],
})
export class RunDetailComponent implements OnInit {
  run?: PrimerRunRecord;
  align?: AlignmentsResponse;
  analysis?: TargetAnalysisResponse;

  // a compact text rendering of GC distribution
  gcSummary: string = '';

  constructor(
    private route: ActivatedRoute,
    private api: PrimersService,
    private seqStore: SequenceSharedService
  ) {}

  ngOnInit(): void {
    const runId = this.route.snapshot.paramMap.get('id');
    if (!runId) return;

    this.api.getRun(runId).subscribe({
      next: (r) => {
        this.run = r;
        const seq = this.seqStore.snapshot.sequence;
        const start = r.targetStart;
        const end = r.targetEnd;
        const target = seq ? seq.substring(start - 1, end) : '';

        if (r.result && target) {
          const req = {
            forward: r.result.forwardPrimer.sequence,
            reverse: r.result.reversePrimer.sequence,
            target: target,
          };
          this.api.alignments(req).subscribe({
            next: (aln) => (this.align = aln),
          });

          const taReq: TargetAnalysisRequest = {
            sequence: target,
            parameters: {
              targetGCMin: 35,
              targetGCMax: 65,
              targetSlidingWindowSize: 20,
              targetHomopolymerMax: 4,
              targetRepeatMax: 5,
            },
          };
          this.api.analyze(taReq).subscribe({
            next: (res) => {
              this.analysis = res;
              this.gcSummary = (res.gcDistribution || []).map((v) => `${v.toFixed(1)}%`).join(' · ');
            },
          });
        }
      },
    });
  }
}
