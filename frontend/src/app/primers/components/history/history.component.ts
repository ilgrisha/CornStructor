// File: frontend/src/app/primers/components/history/history.component.ts
// Version: v0.3.0
/**
 * HistoryComponent (Primers)
 * --------------------------
 * Non-standalone component (declared in PrimersModule) that lists recent primer design runs
 * via PrimersService.listRuns(). Clicking a row navigates to /primers/runs/:id.
 */
import { Component, OnInit } from '@angular/core';
import { PrimerRunRecord } from '../../models/primer-run.model';
import { PrimersService } from '../../services/primers.service';
import { Router } from '@angular/router';

@Component({
  selector: 'app-primers-history',
  templateUrl: './history.component.html',
  styleUrls: ['./history.component.scss'],
})
export class HistoryComponent implements OnInit {
  loading = false;
  error: string | null = null;
  runs: PrimerRunRecord[] = [];

  constructor(private api: PrimersService, private router: Router) {}

  ngOnInit(): void {
    this.fetch();
  }

  /** Refresh the runs list. */
  fetch(): void {
    this.loading = true;
    this.error = null;
    this.api.listRuns().subscribe({
      next: (data) => {
        this.runs = data;
        this.loading = false;
      },
      error: (err) => {
        this.error = err?.error?.detail ?? 'Failed to load runs.';
        this.loading = false;
      },
    });
  }

  /** Navigate to a specific run's report. */
  open(id: string): void {
    this.router.navigate(['/primers/runs', id]);
  }
}
