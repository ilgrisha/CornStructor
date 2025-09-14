// File: frontend/src/app/features/runs-history/runs-history.component.ts
// Version: v0.7.1
/**
 * RunsHistoryComponent
 * Modal listing previous runs with infinite scroll + head polling to
 * prepend newly created/completed runs while the dialog is open.
 *
 * v0.7.1:
 *  - Ensure the deleted run disappears immediately and the list refreshes
 *    from the server after deletion succeeds (no need to close/reopen).
 *    * Optimistic removal from `items` for instant UI feedback.
 *    * On success: call `resetAndLoad()` to re-sync pagination/total.
 *    * On error: restore the removed row.
 *
 * v0.7.0:
 *  - Always refetch first page on modal open + immediate head check.
 */
import {
  AfterViewInit,
  Component,
  ElementRef,
  EventEmitter,
  Input,
  OnInit,
  OnChanges,
  OnDestroy,
  Output,
  SimpleChanges,
  ViewChild,
  signal,
} from '@angular/core';
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
  styleUrls: ['./runs-history.component.css'],
})
export class RunsHistoryComponent implements OnInit, OnChanges, AfterViewInit, OnDestroy {
  @Input() open = false;
  @Output() close = new EventEmitter<void>();
  /** Fired when the user chooses a run (job_id). */
  @Output() loadDesign = new EventEmitter<string>();

  @ViewChild('scrollHost') scrollHost?: ElementRef<HTMLElement>;
  @ViewChild('sentinel') sentinel?: ElementRef<HTMLElement>;

  /** Search text. */
  search = '';

  /** List state. */
  readonly items = signal<RunItem[]>([]);
  readonly loadingFirst = signal(false);
  readonly loadingMore = signal(false);
  readonly hasMore = signal(true);

  /** Paging. */
  private readonly pageSize = 50;
  private offset = 0;
  private total = 0;

  /** IO for infinite scroll. */
  private io?: IntersectionObserver;

  /** Poll timer id. */
  private pollId?: number;

  constructor(
    private runs: RunsService,
    private designs: DesignService,
    private analysis: AnalysisService
  ) {}

  ngOnInit(): void {
    if (this.open) this.resetAndLoad();
  }

  ngOnChanges(changes: SimpleChanges): void {
    if (changes['open']) {
      const nowOpen = !!changes['open'].currentValue;
      const wasOpen = !!changes['open'].previousValue;

      if (nowOpen && !wasOpen) {
        this.resetAndLoad();
        this.startPolling(true /* immediateHeadCheck */);
      } else if (!nowOpen && wasOpen) {
        this.stopPolling();
        this.disconnectIO();
      }
    }
  }

  ngAfterViewInit(): void {
    this.setupIO();
  }

  ngOnDestroy(): void {
    this.stopPolling();
    this.disconnectIO();
  }

  /** Handle search action. */
  onSearch() {
    this.resetAndLoad();
  }

  /** Clear list and fetch first page. */
  private resetAndLoad() {
    this.items.set([]);
    this.offset = 0;
    this.total = 0;
    this.hasMore.set(true);
    this.loadingFirst.set(true);

    const q = this.search.trim().length ? this.search.trim() : null;
    this.runs.list(q, this.pageSize, this.offset).subscribe({
      next: (res) => {
        this.total = res.total ?? 0;
        this.items.set(res.items ?? []);
        this.offset = this.items().length;
        this.hasMore.set(this.offset < this.total);
        this.loadingFirst.set(false);
        queueMicrotask(() => this.observeSentinel());
      },
      error: () => {
        this.loadingFirst.set(false);
        this.hasMore.set(false);
      },
    });
  }

  /** Load the next page and append. */
  private loadMore() {
    if (this.loadingMore() || !this.hasMore()) return;

    this.loadingMore.set(true);
    const q = this.search.trim().length ? this.search.trim() : null;

    this.runs.list(q, this.pageSize, this.offset).subscribe({
      next: (res) => {
        const current = this.items();
        const seen = new Set(current.map((i) => i.job_id));
        const incoming = (res.items ?? []).filter((i) => !seen.has(i.job_id));
        this.items.set([...current, ...incoming]);

        this.offset = this.items().length;
        this.total = res.total ?? this.total;
        this.hasMore.set(this.offset < this.total);
        this.loadingMore.set(false);
      },
      error: () => {
        this.loadingMore.set(false);
      },
    });
  }

  /** Poll the head (first page) to prepend new runs while the dialog is open. */
  private startPolling(immediateHeadCheck = false) {
    this.stopPolling();
    if (immediateHeadCheck) this.checkForNewHead();
    this.pollId = window.setInterval(() => this.checkForNewHead(), 4000);
  }

  private stopPolling() {
    if (this.pollId) {
      window.clearInterval(this.pollId);
      this.pollId = undefined;
    }
  }

  /** Fetch first page; prepend any unseen runs and update statuses of existing ones. */
  private checkForNewHead() {
    if (!this.open || this.loadingFirst()) return;

    const q = this.search.trim().length ? this.search.trim() : null;
    this.runs.list(q, this.pageSize, 0).subscribe({
      next: (res) => {
        const current = this.items();
        const firstPage = res.items ?? [];

        if (!current.length) {
          this.items.set(firstPage);
          this.offset = this.items().length;
          this.total = res.total ?? 0;
          this.hasMore.set(this.offset < this.total);
          return;
        }

        const byId = new Map(current.map((i) => [i.job_id, i]));

        // Update existing
        for (const item of firstPage) {
          const existing = byId.get(item.job_id);
          if (existing) {
            if (
              existing.status !== item.status ||
              existing.report_url !== item.report_url ||
              existing.updated_at !== item.updated_at
            ) {
              byId.set(item.job_id, { ...existing, ...item });
            }
          }
        }

        // New at head
        const existingIds = new Set(byId.keys());
        const fresh = firstPage.filter((i) => !existingIds.has(i.job_id));

        const updatedCurrent = Array.from(byId.values()).sort(
          (a, b) => new Date(b.created_at).getTime() - new Date(a.created_at).getTime()
        );

        this.items.set([...fresh, ...updatedCurrent]);
        this.offset = this.items().length;
        this.total = res.total ?? this.total;
        this.hasMore.set(this.offset < this.total);
      },
      error: () => {
        /* ignore transient errors */
      },
    });
  }

  /** Row actions */
  select(r: RunItem) {
    this.loadDesign.emit(r.job_id);
  }

  onDelete(r: RunItem, ev?: MouseEvent) {
    // prevent row click selection when pressing the button
    if (ev) ev.stopPropagation();

    const ok = window.confirm(`Delete run ${r.job_id}? This will also delete its design and artifacts.`);
    if (!ok) return;

    // 1) Optimistic UI: remove from current list immediately
    const before = this.items();
    const removed = before.find((x) => x.job_id === r.job_id);
    this.items.set(before.filter((x) => x.job_id !== r.job_id));

    // 2) Call API
    this.runs.delete(r.job_id, true).subscribe({
      next: () => {
        // 3) Re-sync from server so totals/pagination are correct
        this.resetAndLoad();
      },
      error: () => {
        // 4) On error, restore the removed row
        if (removed) {
          this.items.set([removed, ...this.items()]);
        }
      },
    });
  }

  openReport(r: RunItem) {
    if (r.report_url) window.open(r.report_url, '_blank', 'noopener');
  }

  trackByJob(_: number, item: RunItem) {
    return item.job_id;
  }

  /** IO setup */
  private setupIO() {
    if (this.io) this.io.disconnect();
    const root = this.scrollHost?.nativeElement || null;
    if (!root) return;

    this.io = new IntersectionObserver(
      (entries) => {
        for (const e of entries) {
          if (e.isIntersecting) this.loadMore();
        }
      },
      { root, rootMargin: '200px', threshold: 0.01 }
    );

    this.observeSentinel();
  }

  private observeSentinel() {
    if (!this.io) return;
    const el = this.sentinel?.nativeElement;
    if (el) {
      this.io.unobserve(el);
      this.io.observe(el);
    }
  }

  private disconnectIO() {
    if (this.io) {
      this.io.disconnect();
      this.io = undefined;
    }
  }
}
