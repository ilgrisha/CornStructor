// File: frontend/src/app/shared/components/go-to-primers-button/go-to-primers-button.component.ts
// Version: v0.2.0
/**
 * GoToPrimersButtonComponent
 * - Shares sequence/start/end to SequenceSharedService
 * - Logs the handoff
 * - Navigates to /primers/design (with robust navigation)
 */
import { Component, Input } from '@angular/core';
import { Router } from '@angular/router';
import { SequenceSharedService } from '../../../shared/services/sequence-shared.service';
import { LogBusService } from '../../../shared/services/log-bus.service';

@Component({
  selector: 'app-go-to-primers-button',
  templateUrl: './go-to-primers-button.component.html',
  styleUrls: ['./go-to-primers-button.component.scss'],
})
export class GoToPrimersButtonComponent {
  @Input() sequence = '';
  @Input() start: number | null = null;  // 1-based inclusive
  @Input() end: number | null = null;    // 1-based inclusive
  @Input() name: string | null = null;

  constructor(
    private router: Router,
    private seqStore: SequenceSharedService,
    private log: LogBusService
  ) {}

  async go(ev?: Event) {
    ev?.preventDefault();
    ev?.stopPropagation();

    if (!this.sequence) {
      this.log.append('Primers: No sequence provided from Construction Tree.');
      return;
    }
    if (this.start == null || this.end == null || this.end < this.start) {
      this.log.append('Primers: Invalid start/end provided from Construction Tree.');
      return;
    }

    // Share state
    this.seqStore.set({
      sequence: this.sequence,
      start: this.start,
      end: this.end,
      name: this.name ?? null,
    });
    this.log.append(
      `Primers: Received "${this.name ?? 'sequence'}" (${this.start}–${this.end}, len ${this.end - this.start + 1}).`
    );

    // Robust navigation attempt
    try {
      const ok = await this.router.navigateByUrl('/primers/design', { replaceUrl: false });
      if (!ok) {
        this.log.append('Primers: navigation to /primers/design failed (router returned false).');
      }
    } catch (e: any) {
      this.log.append(`Primers: navigation error → ${e?.message ?? e}`);
    }
  }
}
