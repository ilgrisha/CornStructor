// File: frontend/src/app/features/runs-history/runs-history.component.ts
// Version: v0.2.0
/**
 * RunsHistoryComponent
 * Modal listing previous runs pulled from /api/runs.
 *
 * - Uses `search: string` bound via [(ngModel)].
 * - On open (at init or later), auto-refreshes the list.
 */
import { Component, EventEmitter, Input, OnChanges, OnInit, Output, SimpleChanges, signal } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { RunsService, RunItem } from '../../core/services/runs.service';

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

  /** Search query bound to the input field. */
  search = '';

  /** Items returned by the API. */
  items = signal<RunItem[]>([]);

  /** Loading state for the table. */
  loading = signal<boolean>(false);

  constructor(private runs: RunsService) {}

  ngOnInit() {
    if (this.open) this.refresh();
  }

  ngOnChanges(changes: SimpleChanges): void {
    // If the modal just opened after initial render, fetch data.
    if (changes['open'] && this.open) {
      this.refresh();
    }
  }

  /** Fetch runs from the API, applying an optional free-text filter. */
  refresh() {
    this.loading.set(true);
    const q = this.search.trim() || null;
    this.runs.list(q, 100, 0).subscribe({
      next: res => {
        this.items.set(res.items);
        this.loading.set(false);
      },
      error: _ => this.loading.set(false),
    });
  }

  /** Open the run's report in a new tab (if available). */
  openReport(r: RunItem) {
    if (r.report_url) window.open(r.report_url, '_blank', 'noopener');
  }
}
