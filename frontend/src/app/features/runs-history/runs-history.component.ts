// File: frontend/src/app/features/runs-history/runs-history.component.ts
// Version: v0.6.0
/**
 * RunsHistoryComponent
 * Modal listing previous runs with infinite scroll + lightweight polling to
 * prepend newly created/completed runs while the dialog is open.
 *
 * v0.6.0:
 *  - Add head polling (every 4s) to detect and prepend newly added runs
 *    without disrupting scroll or requiring user action.
 *  - Keep existing infinite scroll "load more" behavior for older pages.
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

  ngOnChanges(_: SimpleChanges): void {
    if (this.open) {
      if (!this.items().length) this.resetAndLoad();
      this.startPolling();
    } else {
      this.stopPolling();
      this.disconnectIO();
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
        this.total = res.total || 0;
        this.items.set(res.items || []);
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
        const incoming = (res.items || []).filter((i) => !seen.has(i.job_id));
        this.items.set([...current, ...incoming]);

        this.offset = this.items().length;
        this.total = res.total || this.total;
        this.hasMore.set(this.offset < this.total);
        this.loadingMore.set(false);
      },
      error: () => {
        this.loadingMore.set(false);
      },
    });
  }

  /** Poll the head (first page) to prepend new runs while the dialog is open. */
  private startPolling() {
    this.stopPolling();
    this.pollId = window.setInterval(() => this.checkForNewHead(), 4000);
  }

  private stopPolling() {
    if (this.pollId) {
      window.clearInterval(this.pollId);
      this.pollId = undefined;
    }
  }

  /** Fetch first page; prepend any unseen runs (by job_id). */
  private checkForNewHead() {
    if (!this.open || this.loadingFirst()) return;

    const q = this.search.trim().length ? this.search.trim() : null;
    this.runs.list(q, this.pageSize, 0).subscribe({
      next: (res) => {
        const current = this.items();
        if (!current.length) {
          // If list is empty (edge), just adopt the page
          this.items.set(res.items || []);
          this.offset = this.items().length;
          this.total = res.total || 0;
          this.hasMore.set(this.offset < this.total);
          return;
        }

        const existingIds = new Set(current.map((i) => i.job_id));
        const fresh = (res.items || []).filter((i) => !existingIds.has(i.job_id));

        if (fresh.length) {
          // Prepend fresh items, keeping overall order (fresh already newest first)
          this.items.set([...fresh, ...current]);
          // Bump counters
          this.offset = this.items().length;
          this.total = res.total || this.total;
          this.hasMore.set(this.offset < this.total);

          // Optional: keep scroll position stable if near top
          const host = this.scrollHost?.nativeElement;
          if (host && host.scrollTop < 40) {
            // stay at top
          }
        }
      },
      error: () => {
        // ignore transient errors
      },
    });
  }

  /** IntersectionObserver setup to watch the sentinel within the scroll container. */
  private setupIO() {
    if (this.io) this.io.disconnect();
    const root = this.scrollHost?.nativeElement || null;
    if (!root) return;

    this.io = new IntersectionObserver(
      (entries) => {
        for (const e of entries) {
          if (e.isIntersecting) {
            this.loadMore();
          }
        }
      },
      {
        root,
        rootMargin: '200px', // start loading a bit before reaching the end
        threshold: 0.01,
      }
    );

    this.observeSentinel();
  }

  private observeSentinel() {
    if (!this.io) return;
    const el = this.sentinel?.nativeElement;
    if (el) {
      this.io.unobserve(el); // idempotent
      this.io.observe(el);
    }
  }

  private disconnectIO() {
    if (this.io) {
      this.io.disconnect();
      this.io = undefined;
    }
  }

  /** Row actions */
  select(r: RunItem) {
    this.loadDesign.emit(r.job_id);
  }

  onDelete(r: RunItem) {
    const ok = window.confirm(`Delete run ${r.job_id}? This will also delete its design and artifacts.`);
    if (!ok) return;
    // Optimistic UI: remove row immediately
    this.items.set(this.items().filter((x) => x.job_id !== r.job_id));
    this.runs.delete(r.job_id, true).subscribe({
      next: () => {
        // adjust counts
        this.total = Math.max(0, this.total - 1);
        this.hasMore.set(this.offset < this.total);
      },
      error: () => {
        // If delete fails, refresh to resync
        this.resetAndLoad();
      },
    });
  }

  openReport(r: RunItem) {
    if (r.report_url) window.open(r.report_url, '_blank', 'noopener');
  }

  trackByJob(_: number, item: RunItem) {
    return item.job_id;
  }
}
