// File: frontend/src/app/app-routing.module.ts
// Version: v0.1.5
/**
 * AppRoutingModule
 * ----------------
 * Top-level routes for CornStructor:
 *  • "/"            → HomeComponent (main landing: Parameters / Sequence / Construction Tree)
 *  • "/run"         → RunComponent (Construction tree design screen)
 *  • "/history"     → RunsHistoryComponent (open CT run history directly)
 *  • "/parameters"  → ParametersEditorComponent (open CT parameter editor as a page)
 *  • "/primers/**"  → lazy-loaded Primers feature (design, history, run detail)
 *  • "/home"        → redirects to "/"
 *  • "**"           → fallback to "/"
 *
 * Notes:
 *  - HomeComponent, RunComponent, RunsHistoryComponent, ParametersEditorComponent are standalone;
 *    we lazy-load them via `loadComponent`.
 *  - PrimersModule remains a lazy-loaded feature module under /primers.
 */
import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';

const routes: Routes = [
  {
    path: '',
    loadComponent: () =>
      import('./features/home/home.component').then((m) => m.HomeComponent),
  },
  {
    path: 'run',
    loadComponent: () =>
      import('./features/run/run.component').then((m) => m.RunComponent),
  },
  {
    path: 'history',
    loadComponent: () =>
      import('./features/runs-history/runs-history.component').then((m) => m.RunsHistoryComponent),
  },
  {
    path: 'parameters',
    loadComponent: () =>
      import('./features/parameters-editor/parameters-editor.component').then((m) => m.ParametersEditorComponent),
  },
  {
    path: 'primers',
    loadChildren: () =>
      import('./primers/primers.module').then((m) => m.PrimersModule),
  },
  { path: 'home', pathMatch: 'full', redirectTo: '' },
  { path: '**', redirectTo: '' },
];

@NgModule({
  imports: [RouterModule.forRoot(routes, { scrollPositionRestoration: 'enabled' })],
  exports: [RouterModule],
})
export class AppRoutingModule {}
