// File: frontend/src/app/features/run/run.component.ts
// Version: v2.7.0
/**
 * RunComponent
 * - Single shared log box: removes the local log pane and only shows SharedLogViewComponent.
 * - Keeps mirroring Construction Tree logs into LogBusService so the Primers tab sees the same stream.
 * - Keeps the Go to Primers button in the header.
 */
import { Component, signal, effect, AfterViewInit, computed } from '@angular/core';
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
  historyOpen = signal<boolean>(false);
  paramsOpen = signal<boolean>(false);

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

  /** Expose sequence/name/start/end to Primers button */
  sequenceText = computed(() => this.a.sequence() ?? '');
  sequenceName = computed(() => this.deriveSequenceName());
  selectedStart = computed(() => this.deriveStart());
  selectedEnd = computed(() => this.deriveEnd());

  constructor(
    public api: ApiService,
    private a: AnalysisService,
    public tp: TreeParamsService,
    private logBus: LogBusService
  ) {}

  ngAfterViewInit(): void {}

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

  openHistory() { this.historyOpen.set(true); }
  closeHistory() { this.historyOpen.set(false); }

  openParameters() {
    this.paramsOpen.set(true);
    if (!this.tp.globals() || Object.keys(this.tp.globals()).length === 0) {
      this.tp.loadDefaults();
    }
  }
  closeParameters() { this.paramsOpen.set(false); }

  private deriveSequenceName(): string {
    const seq = this.a.sequence() ?? '';
    if (!seq) return 'sequence';
    const anyA: any = this.a as any;
    const header = typeof anyA.header === 'function' ? anyA.header() : (anyA.header ?? null);
    if (header && typeof header === 'string' && header.trim()) return header.trim();
    return `sequence (${seq.length} nt)`;
  }

  private deriveStart(): number {
    const g = this.tp.globals?.();
    const seq = this.a.sequence() ?? '';
    const byGlobals = (g && (g as any).target && Number((g as any).target.start)) || null;
    const start = byGlobals ?? 1;
    return Math.max(1, Math.min(start, Math.max(1, seq.length)));
  }

  private deriveEnd(): number {
    const g = this.tp.globals?.();
    const seq = this.a.sequence() ?? '';
    const byGlobals = (g && (g as any).target && Number((g as any).target.end)) || null;
    const endDefault = seq ? seq.length : 1;
    const end = byGlobals ?? endDefault;
    return Math.max(1, Math.min(end, Math.max(1, seq.length)));
  }
}
