// File: frontend/src/app/shared/services/log-bus.service.ts
// Version: v0.1.0
/**
 * LogBusService
 * -------------
 * Lightweight event bus to share log messages across tabs (Construction Tree & Primers).
 * Append messages during long-running tasks; consumers can render a shared Log box.
 */

import { Injectable } from '@angular/core';
import { BehaviorSubject } from 'rxjs';

@Injectable({ providedIn: 'root' })
export class LogBusService {
  private lines$ = new BehaviorSubject<string[]>([]);

  get stream() {
    return this.lines$.asObservable();
  }

  append(line: string) {
    this.lines$.next([...this.lines$.value, line]);
  }

  appendBlock(block: string[]) {
    this.lines$.next([...this.lines$.value, ...block]);
  }

  clear() {
    this.lines$.next([]);
  }
}
