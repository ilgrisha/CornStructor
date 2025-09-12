/* ============================================================================
 * Path: frontend/src/app/core/services/api.service.ts
 * Version: v2.3.0
 * ============================================================================
 * - startDesign() posts to backend and streams logs via SSE
 * - Signals for logs + result links (supports absolute or relative URLs)
 * - Looks for a final RESULT pointing to /reports/{jobId}/index.html
 * ==========================================================================*/
import { Injectable, signal, WritableSignal } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { AnalysisParams } from './analysis.service';

export interface DesignStartRequest {
  sequence: string;
  params: AnalysisParams;
  toggles?: Record<string, boolean>;
}

@Injectable({ providedIn: 'root' })
export class ApiService {
  private sse?: EventSource;

  // UI signals
  readonly logLines: WritableSignal<string[]> = signal<string[]>([]);
  readonly resultLinks: WritableSignal<string[]> = signal<string[]>([]);
  readonly resultIndexLink: WritableSignal<string | null> = signal<string | null>(null);

  constructor(private http: HttpClient) {}

  startDesign(req: DesignStartRequest) {
    this.dispose();
    this.logLines.set([]);
    this.resultLinks.set([]);
    this.resultIndexLink.set(null);
    this.http.post<{ jobId: string }>('/api/design/start', req)
      .subscribe({
        next: (res) => this.attachSSE(res.jobId),
        error: (err) => this.logLines.update(a => [...a, `ERROR: ${err?.message ?? err}`])
      });
  }

  private dispose() {
    if (this.sse) {
      this.sse.close();
      this.sse = undefined;
    }
  }

  private attachSSE(jobId: string) {
    this.sse = new EventSource(`/api/design/${jobId}/logs`);
    this.sse.onmessage = (ev) => {
      const line = ev.data as string;
      this.logLines.update(a => [...a, line]);

      // Capture RESULT links. Accept absolute (http/https) or relative (/reports/...)
      const m = /^RESULT:\s*(\S+)/.exec(line);
      if (m) {
        const url = m[1];
        this.resultLinks.update(arr => [...arr, url]);
        if (url.endsWith('/index.html')) {
          this.resultIndexLink.set(url);
        }
      }
    };
    this.sse.onerror = () => {
      this.logLines.update(a => [...a, 'Log stream ended']);
      this.dispose();
    };
  }
}
