// File: frontend/src/app/shared/components/go-to-primers-button/go-to-primers-button.component.ts
// Version: v0.1.0
/**
 * GoToPrimersButtonComponent
 * --------------------------
 * A drop-in button that:
 *  1) Saves the current Sequence + target window into SequenceSharedService
 *  2) Appends a message to the shared Log
 *  3) Navigates to /primers/design
 *
 * Usage (in Construction Tree template):
 * <app-go-to-primers-button
 *   [sequence]="sequenceText"
 *   [start]="selectedStart"
 *   [end]="selectedEnd"
 *   [name]="sequenceName">
 * </app-go-to-primers-button>
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
  /**
   * Full raw sequence (no FASTA header). If you have FASTA, strip headers before passing.
   */
  @Input() sequence = '';

  /**
   * 1-based inclusive start coordinate (required to prefill the Primers tab).
   */
  @Input() start: number | null = null;

  /**
   * 1-based inclusive end coordinate (required to prefill the Primers tab).
   */
  @Input() end: number | null = null;

  /**
   * Optional display name for the sequence (e.g., FASTA header or record id).
   */
  @Input() name: string | null = null;

  constructor(
    private router: Router,
    private seqStore: SequenceSharedService,
    private log: LogBusService
  ) {}

  go() {
    if (!this.sequence) {
      this.log.append('Primers: No sequence provided from Construction Tree.');
      return;
    }
    if (this.start == null || this.end == null || this.end < this.start) {
      this.log.append('Primers: Invalid start/end provided from Construction Tree.');
      return;
    }

    // Share sequence state
    this.seqStore.set({
      sequence: this.sequence,
      start: this.start,
      end: this.end,
      name: this.name ?? null,
    });

    // Log the transition
    this.log.append(
      `Primers: Launching design for ${this.name ?? 'sequence'} ` +
      `(${this.start}–${this.end}, length ${this.end - this.start + 1}).`
    );

    // Navigate to the Primers Design tab
    this.router.navigate(['/primers/design']);
  }
}
