// File: frontend/src/app/features/run/run.component.ts
// Version: v2.6.2
/**
 * RunComponent (updated)
 * ----------------------
 * Fix: since this component is standalone, import non-standalone shared components via SharedModule.
 * - Replaces direct imports of GoToPrimersButtonComponent and SharedLogViewComponent with SharedModule.
 * - Removes unused SequenceSharedService injection.
 *
 * Keeps:
 * - Construction Tree workflow
 * - "Go to Primers" button
 * - Log mirroring to shared LogBusService
 */
import { Component, signal, effect, ViewChild, ElementRef, AfterViewInit, computed } from '@angular/core';
import { CommonModule } from '@angular/common';
import { ApiService } from '../../core/services/api.service';
import { AnalysisService } from '../../core/services/analysis.service';
import { RunsHistoryComponent } from '../runs-history/runs-history.component';
import { ParametersEditorComponent } from '../parameters-editor/parameters-editor.component';
import { TreeParamsService } from '../../core/services/tree-params.service';

// 🔗 Shared module exports the non-standalone components we need
import { SharedModule } from '../../shared/shared.module';
import { LogBusService } from '../../shared/services/log-bus.service';

@Component({
  selector: 'app-run',
  standalone: true,
  imports: [CommonModule, RunsHistoryComponent, ParametersEditorComponent, SharedModule],
  templateUrl: './run.component.html',
  styleUrls: ['./run.component.css']
})
export class RunComponent implements AfterViewInit {
  /** Modals */
  historyOpen = signal<boolean>(false);
  paramsOpen = signal<boolean>(false);

  /** Scrollable logs container */
  @ViewChild('logPane') logPane?: ElementRef<HTMLDivElement>;

  /** Mirror Construction Tree logs into shared Log bus */
  private prevLen = 0;
  private mirrorEffect = effect(() => {
    const lines = this.api.logLines();
    if (lines && lines.length > this.prevLen) {
      const newChunk = lines.slice(this.prevLen);
      this.logBus.appendBlock(newChunk);
      this.prevLen = lines.length;
    }
  });

  /** Auto-scroll the local log pane */
  private autoScrollEffect = effect(() => {
    const lines = this.api.logLines();
    const _len = lines.length;
    queueMicrotask(() => {
      const el = this.logPane?.nativeElement;
      if (el) el.scrollTop = el.scrollHeight;
    });
  });

  /** Expose sequence & target window for the Primers button */
  sequenceText = computed(() => this.a.sequence() ?? '');
  sequenceName = computed(() => this.deriveSequenceName());

  /** Selected coordinates (1-based inclusive); fall back to full sequence if unspecified */
  selectedStart = computed(() => this.deriveStart());
  selectedEnd = computed(() => this.deriveEnd());

  constructor(
    public api: ApiService,
    private a: AnalysisService,
    public tp: TreeParamsService,
    private logBus: LogBusService
  ) {}

  ngAfterViewInit(): void {
    setTimeout(() => {
      const el = this.logPane?.nativeElement;
      if (el) el.scrollTop = el.scrollHeight;
    }, 0);
  }

  /** Start a new Construction Tree design run using current inputs + tree params. */
  run() {
    const req = {
      sequence: this.a.sequence(),
      params: {
        analysis: this.a.params(),
        tree: { globals: this.tp.globals(), levels: this.tp.levels() }
      },
      toggles: { homopolymers: true, gc: true, entropy: true, repeats: true }
    };
    this.api.startDesign(req);
  }

  /** History modal */
  openHistory() { this.historyOpen.set(true); }
  closeHistory() { this.historyOpen.set(false); }

  /** Parameters modal */
  openParameters() {
    this.paramsOpen.set(true);
    if (!this.tp.globals() || Object.keys(this.tp.globals()).length === 0) {
      this.tp.loadDefaults();
    }
  }
  closeParameters() { this.paramsOpen.set(false); }

  // ---- Helpers to derive start/end/name from existing services ----

  /** Try to pull a friendly name for the sequence, otherwise fallback. */
  private deriveSequenceName(): string {
    const seq = this.a.sequence() ?? '';
    if (!seq) return 'sequence';
    const anyA: any = this.a as any;
    const header = typeof anyA.header === 'function' ? anyA.header() : (anyA.header ?? null);
    if (header && typeof header === 'string' && header.trim()) return header.trim();
    return `sequence (${seq.length} nt)`;
  }

  /** Try to read a selected start coordinate from TreeParamsService.globals().target.start (1-based). */
  private deriveStart(): number {
    const g = this.tp.globals?.();
    const seq = this.a.sequence() ?? '';
    const byGlobals = (g && (g as any).target && Number((g as any).target.start)) || null;
    const start = byGlobals ?? 1;
    return Math.max(1, Math.min(start, Math.max(1, seq.length)));
  }

  /** Try to read a selected end coordinate from TreeParamsService.globals().target.end (1-based). */
  private deriveEnd(): number {
    const g = this.tp.globals?.();
    const seq = this.a.sequence() ?? '';
    const byGlobals = (g && (g as any).target && Number((g as any).target.end)) || null;
    const endDefault = seq ? seq.length : 1;
    const end = byGlobals ?? endDefault;
    return Math.max(1, Math.min(end, Math.max(1, seq.length)));
  }
}
