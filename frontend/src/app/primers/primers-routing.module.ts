// File: frontend/src/app/primers/primers-routing.module.ts
// Version: v0.3.1
/**
 * Primers feature routes under /primers
 *   /primers/design
 *   /primers/history
 *   /primers/runs/:id
 */
import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';
import { PrimerDesignComponent } from './components/primer-design/primer-design.component';
import { HistoryComponent } from './components/history/history.component';
import { RunDetailComponent } from './components/run-detail/run-detail.component';

const routes: Routes = [
  { path: 'design', component:  PrimerDesignComponent },
  { path: 'history', component: HistoryComponent },
  { path: 'runs/:id', component: RunDetailComponent },
  { path: '', pathMatch: 'full', redirectTo: 'design' },
];

@NgModule({
  imports: [RouterModule.forChild(routes)],
  exports: [RouterModule],
})
export class PrimersRoutingModule {}
