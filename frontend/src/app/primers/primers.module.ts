// File: frontend/src/app/primers/primers.module.ts
// Version: v0.4.0
/**
 * PrimersModule (updated)
 * -----------------------
 * - Declares non-standalone HistoryComponent and others.
 * - Imports the standalone PrimersParametersEditorComponent so it is usable in templates.
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
import { PrimersParametersEditorComponent } from './components/parameters-editor/parameters-editor.component'; // <-- standalone

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
    SharedModule,
    // bring the standalone modal into this module's scope
    PrimersParametersEditorComponent
  ]
})
export class PrimersModule {}
