// File: frontend/src/app/primers/primers-routing.module.ts
// Version: v0.1.0
/**
 * PrimersRoutingModule
 * --------------------
 * Routes for the primers feature:
 *  - /primers/design
 *  - /primers/runs/:id
 *  - /primers/history
 */
import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';
import { PrimerDesignComponent } from './components/primer-design/primer-design.component';
import { RunDetailComponent } from './components/run-detail/run-detail.component';
import { HistoryComponent } from './components/history/history.component';

const routes: Routes = [
  { path: 'design', component: PrimerDesignComponent },
  { path: 'runs/:id', component: RunDetailComponent },
  { path: 'history', component: HistoryComponent },
  { path: '', pathMatch: 'full', redirectTo: 'design' }
];

@NgModule({
  imports: [RouterModule.forChild(routes)],
  exports: [RouterModule],
})
export class PrimersRoutingModule {}
