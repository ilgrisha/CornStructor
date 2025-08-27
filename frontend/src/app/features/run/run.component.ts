// File: frontend/src/app/features/run/run.component.ts
// Version: v2.3.0
/**
 * RunComponent
 * - Launches the design job and streams logs.
 * - Opens the History modal.
 * - Auto-scrolls the log view as new lines arrive.
 */
import { Component, signal, effect, ViewChild, ElementRef, AfterViewInit } from '@angular/core';
import { CommonModule } from '@angular/common';
import { ApiService } from '../../core/services/api.service';
import { AnalysisService } from '../../core/services/analysis.service';
import { RunsHistoryComponent } from '../runs-history/runs-history.component';

@Component({
  selector: 'app-run',
  standalone: true,
  imports: [CommonModule, RunsHistoryComponent],
  templateUrl: './run.component.html',
  styleUrls: ['./run.component.css']
})
export class RunComponent implements AfterViewInit {
  /** Controls History modal visibility. */
  historyOpen = signal<boolean>(false);

  /** Reference to the scrollable logs container. */
  @ViewChild('logPane') logPane?: ElementRef<HTMLDivElement>;

  /**
   * Auto-scroll effect: whenever log lines change, scroll the log box
   * to the bottom to reveal the latest output.
   */
  private autoScrollEffect = effect(() => {
    // Access the signal to create the reactive dependency.
    const lines = this.api.logLines();
    const _len = lines.length; // use length so effect runs on append

    // Defer until after view updates.
    queueMicrotask(() => {
      const el = this.logPane?.nativeElement;
      if (el) {
        el.scrollTop = el.scrollHeight;
      }
    });
  });

  constructor(public api: ApiService, private a: AnalysisService) {}

  ngAfterViewInit(): void {
    // Ensure the first paint is scrolled to bottom too (usually empty, but harmless).
    setTimeout(() => {
      const el = this.logPane?.nativeElement;
      if (el) el.scrollTop = el.scrollHeight;
    }, 0);
  }

  /** Start a new design run using current inputs. */
  run() {
    const req = {
      sequence: this.a.sequence(),
      params: this.a.params(),
      toggles: { homopolymers: true, gc: true, entropy: true, repeats: true }
    };
    this.api.startDesign(req);
  }

  /** Open history modal. */
  openHistory() { this.historyOpen.set(true); }

  /** Close history modal. */
  closeHistory() { this.historyOpen.set(false); }
}
