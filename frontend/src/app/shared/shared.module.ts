// File: frontend/src/app/shared/shared.module.ts
// Version: v0.2.0
/**
 * SharedModule
 * ------------
 * Declares and exports shared UI components/services for reuse across features.
 * - SharedLogViewComponent
 * - GoToPrimersButtonComponent  ← NEW
 */
import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';
import { SharedLogViewComponent } from './components/shared-log-view/shared-log-view.component';
import { GoToPrimersButtonComponent } from './components/go-to-primers-button/go-to-primers-button.component';
import { RouterModule } from '@angular/router';

@NgModule({
  declarations: [
    SharedLogViewComponent,
    GoToPrimersButtonComponent
  ],
  imports: [
    CommonModule,
    RouterModule
  ],
  exports: [
    SharedLogViewComponent,
    GoToPrimersButtonComponent
  ],
})
export class SharedModule {}
