// File: frontend/src/app/features/runs-history/runs-history.component.ts
// Version: v0.3.0
/**
 * RunsHistoryComponent
 * Modal listing previous runs pulled from /api/runs.
 *
 * New in v0.3.0:
 * - Emits (loadDesign) with the selected run's job_id.
 * - Convenience method loadDesignAndApply() that fetches the Design and updates
 *   the global AnalysisService.sequence so the Sequence Box reflects it.
 */
import { Component, EventEmitter, Input, OnChanges, OnInit, Output, SimpleChanges, signal } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { RunsService, RunItem } from '../../core/services/runs.service';
import { DesignService } from '../../core/services/design.service';
import { AnalysisService } from '../../core/services/analysis.service';

@Component({
  selector: 'app-runs-history',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './runs-history.component.html',
  styleUrls: ['./runs-history.component.css']
})
export class RunsHistoryComponent implements OnInit, OnChanges {
  @Input() open = false;
  @Output() close = new EventEmitter<void>();
  /** Fired when the user chooses a run (job_id). */
  @Output() loadDesign = new EventEmitter<string>();

  search = '';
  loading = signal(false);
  items = signal<RunItem[]>([]);

  constructor(
    private runs: RunsService,
    private designs: DesignService,
    private analysis: AnalysisService
  ) {}

  ngOnInit(): void {
    if (this.open) this.refresh();
  }

  ngOnChanges(_: SimpleChanges): void {
    if (this.open) this.refresh();
  }

  refresh() {
    this.loading.set(true);
    const q = this.search.trim().length ? this.search.trim() : null;
    this.runs.list(q, 100, 0).subscribe({
      next: res => {
        this.items.set(res.items);
        this.loading.set(false);
      },
      error: _ => this.loading.set(false),
    });
  }

  /** Emit selection to parent. */
  select(r: RunItem) {
    this.loadDesign.emit(r.job_id);
  }

  /** Fetch the Design for a run and apply its sequence to the UI immediately. */
  loadDesignAndApply(r: RunItem) {
    this.designs.getByRun(r.job_id).subscribe({
      next: d => {
        this.analysis.sequence.set(d.sequence);
        this.close.emit();
      },
      error: _ => {
        // TODO: surface a toast/snackbar if desired
      }
    });
  }

  /** Open the run's report in a new tab (if available). */
  openReport(r: RunItem) {
    if (r.report_url) window.open(r.report_url, '_blank', 'noopener');
  }
}
