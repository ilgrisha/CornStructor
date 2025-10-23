// File: frontend/src/app/shared/components/shared-log-view/shared-log-view.component.ts
// Version: v0.1.0
/**
 * SharedLogViewComponent
 * ----------------------
 * Renders the shared log stream so multiple tabs show the same log output.
 */
import { Component, OnDestroy, OnInit } from '@angular/core';
import { LogBusService } from '../../services/log-bus.service';
import { Subscription } from 'rxjs';

@Component({
  selector: 'app-shared-log-view',
  templateUrl: './shared-log-view.component.html',
  styleUrls: ['./shared-log-view.component.scss']
})
export class SharedLogViewComponent implements OnInit, OnDestroy {
  lines: string[] = [];
  private sub?: Subscription;

  constructor(private log: LogBusService) {}

  ngOnInit(): void {
    this.sub = this.log.stream.subscribe(ls => {
      this.lines = ls;
      setTimeout(() => this.scrollToBottom(), 0);
    });
  }

  ngOnDestroy(): void {
    this.sub?.unsubscribe();
  }

  private scrollToBottom() {
    const el = document.querySelector('.shared-log') as HTMLElement | null;
    if (el) el.scrollTop = el.scrollHeight;
  }
}
