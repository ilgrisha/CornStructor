// File: frontend/src/app/core/services/runs.service.ts
// Version: v0.1.0
/**
 * RunsService — talks to /api/runs for history.
 */
import { Injectable, signal, WritableSignal } from '@angular/core';
import { HttpClient } from '@angular/common/http';

export interface RunItem {
  job_id: string;
  status: 'running' | 'completed' | 'failed';
  created_at: string;
  updated_at: string;
  report_url?: string | null;
  exit_code?: number | null;
  sequence_len?: number | null;
  note?: string | null;
}

export interface RunListResponse {
  total: number;
  items: RunItem[];
}

@Injectable({ providedIn: 'root' })
export class RunsService {
  constructor(private http: HttpClient) {}

  list(q: string | null = null, limit = 50, offset = 0) {
    const params: any = { limit, offset };
    if (q) params.q = q;
    return this.http.get<RunListResponse>('/api/runs', { params });
  }
}
