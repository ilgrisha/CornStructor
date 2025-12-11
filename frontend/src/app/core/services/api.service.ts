/* ============================================================================
 * Path: frontend/src/app/core/services/api.service.ts
 * Version: v2.4.1
 * ============================================================================
 * - startDesign() posts to backend and streams logs via SSE
 * - Signals for logs + result links (supports absolute or relative URLs)
 * - Looks for a final RESULT pointing to /reports/{jobId}/index.html
 * - NOTE: `params` type is `any` to allow { analysis, tree }
 * ==========================================================================*/
import { Injectable, signal, WritableSignal } from '@angular/core';
import { HttpClient } from '@angular/common/http';

export interface DesignStartRequest {
  sequence: string;
  params: any;
  toggles?: Record<string, boolean>;
  note?: string | null;
}

interface CompletionInfo {
  jobId: string;
  exitCode: number;
  timestamp: number;
}

@Injectable({ providedIn: 'root' })
export class ApiService {
  private sse?: EventSource;

  readonly logLines: WritableSignal<string[]> = signal<string[]>([]);
  readonly resultLinks: WritableSignal<string[]> = signal<string[]>([]);
  readonly resultIndexLink: WritableSignal<string | null> = signal<string | null>(null);
  readonly currentJobId: WritableSignal<string | null> = signal<string | null>(null);
  readonly lastCompletion: WritableSignal<CompletionInfo | null> = signal<CompletionInfo | null>(null);

  constructor(private http: HttpClient) {}

  startDesign(req: DesignStartRequest) {
    this.dispose();
    this.logLines.set([]);
    this.resultLinks.set([]);
    this.resultIndexLink.set(null);
    this.currentJobId.set(null);
    this.lastCompletion.set(null);

    this.http.post<{ jobId: string }>('/api/design/start', req)
      .subscribe({
        next: (res) => {
          this.currentJobId.set(res.jobId);
          this.attachSSE(res.jobId);
        },
        error: (err) => this.logLines.update(a => [...a, `ERROR: ${err?.message ?? err}`])
      });
  }

  private dispose() {
    if (this.sse) {
      this.sse.close();
      this.sse = undefined;
    }
    this.currentJobId.set(null);
  }

  private attachSSE(jobId: string) {
    this.sse = new EventSource(`/api/design/${jobId}/logs`);
    this.sse.onmessage = (ev) => {
      const line = ev.data as string;
      this.logLines.update(a => [...a, line]);

      const m = /^RESULT:\s*(\S+)/.exec(line);
      if (m) {
        const url = m[1];
        this.resultLinks.update(arr => [...arr, url]);
        if (url.endsWith('/index.html')) this.resultIndexLink.set(url);
      }
      const exit = /^EXIT:\s*(-?\d+)/.exec(line);
      if (exit) {
        const code = Number(exit[1]);
        this.lastCompletion.set({ jobId, exitCode: code, timestamp: Date.now() });
        if (code !== 0) this.resultIndexLink.set(null);
        this.currentJobId.set(null);
      }
    };
    this.sse.onerror = () => {
      this.logLines.update(a => [...a, 'Log stream ended']);
      this.dispose();
    };
  }
}
