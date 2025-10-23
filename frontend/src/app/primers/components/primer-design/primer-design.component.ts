// File: frontend/src/app/primers/components/primer-design/primer-design.component.ts
// Version: v0.2.0
/**
 * PrimerDesignComponent (updated)
 * -------------------------------
 * - Parameters button now opens a modal editor (loads via GET, saves via PUT).
 * - Start Design uses the server-stored parameters (do not include params in POST),
 *   keeping parity with Construction Tree "Parameters" UX.
 */
import { Component, OnDestroy, OnInit } from '@angular/core';
import { FormBuilder, FormGroup, Validators } from '@angular/forms';
import { PrimersService } from '../../services/primers.service';
import { SequenceSharedService } from '../../../shared/services/sequence-shared.service';
import { LogBusService } from '../../../shared/services/log-bus.service';
import { Router } from '@angular/router';
import { PrimerDesignRequest, PrimerDesignResponse } from '../../models/primer-params.model';
import { Subscription } from 'rxjs';

@Component({
  selector: 'app-primer-design',
  templateUrl: './primer-design.component.html',
  styleUrls: ['./primer-design.component.scss'],
})
export class PrimerDesignComponent implements OnInit, OnDestroy {
  sequence = '';
  seqName: string | null | undefined = null;

  form!: FormGroup;
  designing = false;
  lastRunId: string | null = null;

  // NEW: parameters modal state
  paramsOpen = false;

  private sub?: Subscription;

  constructor(
    private fb: FormBuilder,
    private api: PrimersService,
    private seqStore: SequenceSharedService,
    private log: LogBusService,
    private router: Router
  ) {}

  ngOnInit(): void {
    this.sub = this.seqStore.value$.subscribe((s) => {
      this.sequence = s.sequence || '';
      this.seqName = s.name ?? null;
      if (s.start != null && s.end != null) {
        this.form.patchValue({ start: s.start, end: s.end }, { emitEvent: false });
      }
    });

    this.form = this.fb.group({
      start: [null, [Validators.required, Validators.min(1)]],
      end: [null, [Validators.required, Validators.min(1)]],
    });
  }

  ngOnDestroy(): void {
    this.sub?.unsubscribe();
  }

  openHistory(): void {
    this.router.navigate(['/primers/history']);
  }

  openParameters(): void {
    this.paramsOpen = true;
  }
  closeParameters(): void {
    this.paramsOpen = false;
  }

  startDesign(): void {
    if (!this.sequence || this.form.invalid) {
      this.log.append('Primer Design: Please provide sequence and valid start/end positions.');
      return;
    }
    const start = Number(this.form.value.start);
    const end = Number(this.form.value.end);
    if (end < start || end > this.sequence.length) {
      this.log.append('Primer Design: Invalid coordinates — ensure 1 ≤ start ≤ end ≤ sequence length.');
      return;
    }

    const body: PrimerDesignRequest = {
      sequence: this.sequence,
      seqTargetStartPosition: start,
      seqTargetEndPosition: end,
      // parameters: omitted → server uses stored parameters (GET/PUT managed by modal)
    } as any;

    this.designing = true;
    this.lastRunId = null;

    this.log.append('Primer Design: Submitting design request using current server parameters...');
    this.api.design(body).subscribe({
      next: (res: PrimerDesignResponse) => {
        this.log.append('Primer Design: Design completed successfully.');
        this.log.append(`Forward primer: ${res.forwardPrimer.sequence} (Tm ${res.forwardPrimer.tm.toFixed(1)}°C, GC ${res.forwardPrimer.gc.toFixed(1)}%)`);
        this.log.append(`Reverse primer: ${res.reversePrimer.sequence} (Tm ${res.reversePrimer.tm.toFixed(1)}°C, GC ${res.reversePrimer.gc.toFixed(1)}%)`);
        if (res.warnings?.length) {
          this.log.append('Warnings:');
          res.warnings.forEach(w => this.log.append(` - ${w}`));
        }
        this.lastRunId = res.runId;
        this.designing = false;
      },
      error: (err) => {
        this.log.append(`Primer Design: Failed — ${err?.error?.detail ?? 'Unknown error'}`);
        this.designing = false;
      }
    });
  }

  openReport(): void {
    if (this.lastRunId) {
      this.router.navigate(['/primers/runs', this.lastRunId]);
    }
  }
}
