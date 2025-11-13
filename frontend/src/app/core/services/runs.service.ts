// File: frontend/src/app/core/services/runs.service.ts
// Version: v0.4.0
/**
 * RunsService: list and manage CornStructor runs.
 *
 * v0.4.0:
 *  - Add cache-busting param `t` to avoid any proxy/browser caching on /api/runs.
 *  - No API shape changes.
 */
import { Injectable } from '@angular/core';
import { HttpClient, HttpParams } from '@angular/common/http';
import { environment } from '../../../environments/environment';
import { Observable } from 'rxjs';

export interface RunItem {
  job_id: string;
  status: 'running' | 'completed' | 'failed' | string;
  sequence_len?: number | null;
  report_url?: string | null;
  created_at: string;
  updated_at: string;
}

export interface RunListResponse {
  items: RunItem[];
  total: number;
}

@Injectable({ providedIn: 'root' })
export class RunsService {
  private api = environment.apiBase || '/api';

  constructor(private http: HttpClient) {}

  /** List runs with pagination; use `offset` to fetch the next page (infinite scroll). */
  list(q: string | null, limit = 50, offset = 0): Observable<RunListResponse> {
    let params = new HttpParams()
      .set('limit', limit)
      .set('offset', offset)
      // cache buster
      .set('t', String(Date.now()));
    if (q) params = params.set('q', q);
    return this.http.get<RunListResponse>(`${this.api}/runs`, { params });
  }

  /** Delete a run and its linked design; optionally keep artifacts by passing deleteReports=false. */
  delete(jobId: string, deleteReports = true): Observable<void> {
    const params = new HttpParams().set('delete_reports', String(deleteReports));
    return this.http.delete<void>(`${this.api}/runs/${jobId}`, { params });
  }
}
