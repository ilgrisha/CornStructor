// File: frontend/src/app/primers/primers.module.ts
// Version: v0.3.0
/**
 * PrimersModule (updated)
 * -----------------------
 * - Registers the new PrimersParametersEditorComponent (standalone; used in template).
 */
import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';
import { PrimersRoutingModule } from './primers-routing.module';
import { ReactiveFormsModule, FormsModule } from '@angular/forms';
import { HttpClientModule } from '@angular/common/http';

import { PrimerDesignComponent } from './components/primer-design/primer-design.component';
import { RunDetailComponent } from './components/run-detail/run-detail.component';
import { HistoryComponent } from './components/history/history.component';
import { SharedModule } from '../shared/shared.module';

// Standalone modal imported directly in template; no need to declare.
// import { PrimersParametersEditorComponent } from './components/parameters-editor/parameters-editor.component';

@NgModule({
  declarations: [
    PrimerDesignComponent,
    RunDetailComponent,
    HistoryComponent
  ],
  imports: [
    CommonModule,
    FormsModule,
    ReactiveFormsModule,
    HttpClientModule,
    PrimersRoutingModule,
    SharedModule
  ]
})
export class PrimersModule {}
