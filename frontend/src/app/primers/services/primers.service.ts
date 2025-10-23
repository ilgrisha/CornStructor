// File: frontend/src/app/primers/services/primers.service.ts
// Version: v0.1.3
/**
 * PrimersService
 * --------------
 * Readable, intuitive API client for primer design, analysis, alignments, run history,
 * and editable server-side primer design parameters.
 */
import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import {
  PrimerDesignRequest,
  PrimerDesignResponse,
  TargetAnalysisRequest,
  TargetAnalysisResponse,
  AlignmentsRequest,
  AlignmentsResponse,
  PrimerDesignParameters,
} from '../models/primer-params.model';
import { Observable } from 'rxjs';
import { PrimerRunRecord } from '../models/primer-run.model';
import { environment } from '../../../environments/environment';

@Injectable({ providedIn: 'root' })
export class PrimersService {
  private readonly base: string = (environment as any).apiBase ?? '';

  constructor(private http: HttpClient) {}

  /** Safely join base URL and path without double slashes. */
  private url(path: string): string {
    if (!this.base) return path;
    return `${this.base.replace(/\/+$/, '')}/${path.replace(/^\/+/, '')}`;
  }

  /** Start a primer design run for the given sequence and target window. */
  design(body: { } & PrimerDesignRequest): Observable<PrimerDesignResponse> {
    return this.http.post<PrimerDesignResponse>(this.url('/api/v1/primers/design'), body);
  }

  /** Run target-region analysis (GC windows, repeats, structure flags). */
  analyze(body: TargetAnalysisRequest): Observable<TargetAnalysisResponse> {
    return this.http.post<TargetAnalysisResponse>(this.url('/api/v1/primers/target-analysis'), body);
  }

  /** Compute ASCII alignments for forward/reverse vs. target, self/cross-dimers. */
  alignments(body: AlignmentsRequest): Observable<AlignmentsResponse> {
    return this.http.post<AlignmentsResponse>(this.url('/api/v1/primers/alignments'), body);
  }

  /** List recent primer design runs (for History view). */
  listRuns(): Observable<PrimerRunRecord[]> {
    return this.http.get<PrimerRunRecord[]>(this.url('/api/v1/primers/runs'));
  }

  /** Fetch a single run (for Report view). */
  getRun(id: string): Observable<PrimerRunRecord> {
    return this.http.get<PrimerRunRecord>(this.url(`/api/v1/primers/runs/${encodeURIComponent(id)}`));
  }

  /** Load the current editable primer design parameters from server. */
  getParameters(): Observable<PrimerDesignParameters> {
    return this.http.get<PrimerDesignParameters>(this.url('/api/v1/primers/parameters'));
  }

  /** Persist new primer design parameters on the server. */
  updateParameters(params: PrimerDesignParameters): Observable<PrimerDesignParameters> {
    return this.http.put<PrimerDesignParameters>(this.url('/api/v1/primers/parameters'), params);
  }
}
