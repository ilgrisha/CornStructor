/* ============================================================================
 * Path: frontend/src/app/core/services/api.service.ts
 * Version: v2.2.0
 * ============================================================================
 * - startDesign() posts to backend and streams logs via SSE
 * - Signals for logs + result link
 * ==========================================================================*/
import { Injectable, signal, WritableSignal, effect } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { AnalysisParams, FeatureKind } from './analysis.service';

export interface DesignStartRequest {
  sequence: string;
  params: AnalysisParams;
  toggles?: Record<string, boolean>;
}

@Injectable({ providedIn: 'root' })
export class ApiService {
  constructor(private http: HttpClient) {}

  logLines: WritableSignal<string[]> = signal<string[]>([]);
  resultLink: WritableSignal<string | null> = signal(null);
  private sse?: EventSource;

  clear() {
    this.logLines.set([]);
    this.resultLink.set(null);
    if (this.sse) { this.sse.close(); this.sse = undefined; }
  }

  startDesign(req: DesignStartRequest) {
    this.clear();
    return this.http.post<{ jobId: string }>('/api/design/start', req).subscribe({
      next: ({ jobId }) => this.attachSSE(jobId),
      error: () => this.logLines.update(a => [...a, 'Error starting design']),
    });
  }

  private attachSSE(jobId: string) {
    this.sse = new EventSource(`/api/design/${jobId}/logs`);
    this.sse.onmessage = (ev) => {
      this.logLines.update(a => [...a, ev.data]);
      // naïve: detect final line with URL
      const m = /RESULT:\s*(https?:\/\/\S+)/.exec(ev.data);
      if (m) this.resultLink.set(m[1]);
    };
    this.sse.onerror = () => {
      this.logLines.update(a => [...a, 'Log stream ended']);
      this.sse?.close();
    };
  }
}
