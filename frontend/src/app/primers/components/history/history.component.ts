// File: frontend/src/app/primers/components/history/history.component.ts
// Version: v0.1.0
/**
 * HistoryComponent
 * ----------------
 * Displays table of prior primer design runs. Clicking a row opens the Report.
 */

import { Component, OnInit } from '@angular/core';
import { PrimersService } from '../../services/primers.service';
import { PrimerRunRecord } from '../../models/primer-run.model';
import { Router } from '@angular/router';

@Component({
  selector: 'app-primers-history',
  templateUrl: './history.component.html',
  styleUrls: ['./history.component.scss'],
})
export class HistoryComponent implements OnInit {
  runs: PrimerRunRecord[] = [];
  loading = false;
  error: string | null = null;

  constructor(private api: PrimersService, private router: Router) {}

  ngOnInit(): void {
    this.loading = true;
    this.api.listRuns().subscribe({
      next: (data) => {
        this.runs = data;
        this.loading = false;
      },
      error: (err) => {
        this.error = err?.error?.detail ?? 'Failed to load history.';
        this.loading = false;
      },
    });
  }

  openRun(run: PrimerRunRecord) {
    this.router.navigate(['/primers/runs', run.id]);
  }
}
