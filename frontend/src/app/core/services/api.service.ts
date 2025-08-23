/* ============================================================================
 * Path: frontend/src/app/core/services/api.service.ts
 * Version: v2.3.0
 * ============================================================================
 * Frontend API glue.
 * - Adds a result link signal + resultLink()/setResultLink()
 * - Clears link & logs at startDesign()
 * - Sets result link opportunistically in getResult()
 * ==========================================================================*/

import { Injectable, WritableSignal, signal } from '@angular/core';
import type { AnalysisParams, Toggles } from './analysis.service';

export interface DesignStartRequest {
  sequence: string;
  params: AnalysisParams;
  toggles: Toggles;
}

export interface DesignStartResponse {
  jobId: string;
}

@Injectable({ providedIn: 'root' })
export class ApiService {
  baseUrl = '/api';

  /** Live log lines from SSE or local messages */
  logLines: WritableSignal<string[]> = signal<string[]>([]);

  /** Currently running job (if any) */
  currentJobId: WritableSignal<string | null> = signal<string | null>(null);

  /** Result page URL (set after job completes) */
  private _resultUrl: WritableSignal<string | null> = signal<string | null>(null);

  /** Template-friendly getter used in run.component.html */
  resultLink(): string | null {
    return this._resultUrl();
  }

  /** Allows components to override/clear the result link */
  setResultLink(url: string | null) {
    this._resultUrl.set(url);
  }

  addLogLine(line: string) {
    this.logLines.update(arr => [...arr, line]);
  }

  async startDesign(req: DesignStartRequest): Promise<DesignStartResponse> {
    // reset UI state for a fresh run
    this.setResultLink(null);
    this.logLines.set([]);

    const res = await fetch(`${this.baseUrl}/design/start`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(req),
    });
    if (!res.ok) {
      this.addLogLine(`Start failed: ${res.status} ${res.statusText}`);
      throw new Error(`startDesign failed: ${res.status}`);
    }
    const data = (await res.json()) as DesignStartResponse;
    this.currentJobId.set(data.jobId);
    this.addLogLine(`Started job ${data.jobId}`);
    return data;
  }

  streamLogs(jobId: string): EventSource {
    const es = new EventSource(`${this.baseUrl}/design/${jobId}/logs`);
    es.onmessage = (ev) => this.addLogLine(ev.data ?? '');
    es.onerror = () => this.addLogLine('⚠️ log stream disconnected');
    return es;
  }

  async getResult(jobId: string): Promise<any> {
    const res = await fetch(`${this.baseUrl}/design/${jobId}/result`);
    if (!res.ok) throw new Error(`result fetch failed: ${res.status}`);
    const data = await res.json();

    // Try a few common fields the backend might provide
    const candidate =
      data?.indexUrl ??
      data?.outputPage ??
      data?.output?.indexUrl ??
      data?.links?.index ??
      null;

    if (candidate) this.setResultLink(candidate);
    return data;
  }
}
