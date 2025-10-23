// File: frontend/src/app/app-routing.module.ts
// Version: v0.1.1
/**
 * AppRoutingModule (updated)
 * -------------------------
 * Adds lazy route for primers feature under /primers.
 * NOTE: If your project already has this route, keep a single definition.
 */
import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';

const routes: Routes = [
  // ... your existing routes
  {
    path: 'primers',
    loadChildren: () => import('./primers/primers.module').then(m => m.PrimersModule)
  },
  { path: '', pathMatch: 'full', redirectTo: 'primers' } // adjust to your desired home
];

@NgModule({
  imports: [RouterModule.forRoot(routes, { scrollPositionRestoration: 'enabled' })],
  exports: [RouterModule]
})
export class AppRoutingModule {}
