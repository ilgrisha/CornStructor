// File: frontend/src/app/core/services/design.service.ts
// Version: v0.1.0
/**
 * DesignService: wraps the CornStructor API.
 * - start(sequence) -> Observable<{ jobId: string }>
 * - streamLogs(jobId) -> Observable<string> (SSE lines)
 */
import { Injectable, NgZone } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { environment } from '../../../environments/environment';
import { Observable } from 'rxjs';
import { map } from 'rxjs/operators';

interface StartResponse {
  job_id?: string;  // snake_case (backend)
  jobId?: string;   // camelCase (compat)
}

@Injectable({ providedIn: 'root' })
export class DesignService {
  private api = environment.apiBase;
  private sse = environment.sseBase;

  constructor(private http: HttpClient, private zone: NgZone) {}

  /** Kick off a design run and normalize the id to camelCase. */
  start(sequence: string, params?: any, toggles?: Record<string, boolean>): Observable<{ jobId: string }> {
    const body = { sequence, params, toggles };
    return this.http.post<StartResponse>(`${this.api}/design/start`, body).pipe(
      map(res => {
        const jobId = res.jobId ?? res.job_id;
        if (!jobId) {
          throw new Error('API did not return a jobId');
        }
        return { jobId };
      })
    );
  }

  /**
   * Open an SSE stream for a given jobId.
   * Call `subscription.unsubscribe()` to close the stream.
   */
  streamLogs(jobId: string): Observable<string> {
    if (!jobId) {
      throw new Error('jobId is required for log streaming');
    }
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
}
