// File: frontend/src/app/core/services/design.service.ts
// Version: v0.2.0
/**
 * DesignService: wraps the CornStructor API.
 * - start(sequence, params, toggles) -> Observable<{ jobId: string }>
 * - streamLogs(jobId) -> Observable<string> (SSE lines)
 * - getByRun(jobId) -> Observable<DesignByRun>
 */
import { Injectable, NgZone } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { environment } from '../../../environments/environment';
import { Observable } from 'rxjs';
import { map } from 'rxjs/operators';
import { DesignByRun } from '../models/design';

interface StartResponse {
  job_id?: string;  // snake_case (backend)
  jobId?: string;   // camelCase (compat)
}

@Injectable({ providedIn: 'root' })
export class DesignService {
  private api = environment.apiBase || '/api';
  private sse = environment.sseBase || '/api';

  constructor(private http: HttpClient, private zone: NgZone) {}

  /** Kick off a design run and normalize the id to camelCase. */
  start(sequence: string, params?: any, toggles?: Record<string, boolean>): Observable<{ jobId: string }> {
    const body = { sequence, params, toggles };
    return this.http.post<StartResponse>(`${this.api}/design/start`, body).pipe(
      map(res => ({ jobId: (res.jobId || res.job_id)! }))
    );
  }

  /** Stream text/event-stream logs for a job. */
  streamLogs(jobId: string): Observable<string> {
    const url = `${this.sse}/design/${jobId}/logs`;
    return new Observable<string>(subscriber => {
      const es = new EventSource(url);
      es.onmessage = (ev: MessageEvent) => {
        // Ensure change detection runs inside Angular zone
        this.zone.run(() => subscriber.next(ev.data));
      };
      es.onerror = () => {
        this.zone.run(() => subscriber.error(new Error('SSE connection error')));
        es.close();
      };
      return () => es.close();
    });
  }

  /** Fetch the persisted Design linked to a Run. */
  getByRun(jobId: string) {
    return this.http.get<DesignByRun>(`${this.api}/designs/by-run/${jobId}`);
  }
}
