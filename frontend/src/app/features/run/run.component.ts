// File: frontend/src/app/features/run/run.component.ts
// Version: v2.5.1
/**
 * RunComponent
 * - Launches the design job and streams logs.
 * - Shows a single "Open report" button when an /index.html RESULT appears.
 * - Opens the History and Parameters modals.
 * - Auto-scrolls the log view as new lines arrive.
 */
import { Component, signal, effect, ViewChild, ElementRef, AfterViewInit } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { ApiService } from '../../core/services/api.service';
import { AnalysisService } from '../../core/services/analysis.service';
import { RunsHistoryComponent } from '../runs-history/runs-history.component';
import { ParametersEditorComponent } from '../parameters-editor/parameters-editor.component';
import { TreeParamsService } from '../../core/services/tree-params.service';

@Component({
  selector: 'app-run',
  standalone: true,
  imports: [CommonModule, FormsModule, RunsHistoryComponent, ParametersEditorComponent],
  templateUrl: './run.component.html',
  styleUrls: ['./run.component.css']
})
export class RunComponent implements AfterViewInit {
  /** Modals */
  historyOpen = signal<boolean>(false);
  paramsOpen = signal<boolean>(false);
  runDescription = '';

  /** Scrollable logs container */
  @ViewChild('logPane') logPane?: ElementRef<HTMLDivElement>;

  /** Auto-scroll the log view when new lines arrive. */
  private autoScrollEffect = effect(() => {
    const lines = this.api.logLines();
    const _len = lines.length;
    queueMicrotask(() => {
      const el = this.logPane?.nativeElement;
      if (el) el.scrollTop = el.scrollHeight;
    });
  });

  constructor(public api: ApiService, private a: AnalysisService, public tp: TreeParamsService) {}

  ngAfterViewInit(): void {
    setTimeout(() => {
      const el = this.logPane?.nativeElement;
      if (el) el.scrollTop = el.scrollHeight;
    }, 0);
  }

  /** Start a new design run using current inputs + tree params (Globals & Levels). */
  run() {
    const note = this.runDescription.trim();
    const req = {
      sequence: this.a.sequence(),
      params: {
        analysis: this.a.params(),
        tree: { globals: this.tp.globals(), levels: this.tp.levels() }
      },
      toggles: { homopolymers: true, gc: true, entropy: true, repeats: true },
      note: note.length ? note : undefined
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
}
