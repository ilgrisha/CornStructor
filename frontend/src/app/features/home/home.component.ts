/* ============================================================================
 * Path: frontend/src/app/features/home/home.component.ts
 * Version: v2.0.1
 * ============================================================================
 * FIX: Use selector 'app-root' so it matches the default element in index.html.
 * ==========================================================================*/
import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { SequenceBoxComponent } from '../sequence-box/sequence-box.component';
import { ParamsComponent } from '../params/params.component';
import { ResultsComponent } from '../results/results.component';
import { RunComponent } from '../run/run.component';

@Component({
  selector: 'app-root', // <-- match <app-root> in index.html
  standalone: true,
  imports: [CommonModule, SequenceBoxComponent, ParamsComponent, ResultsComponent, RunComponent],
  templateUrl: './home.component.html',
  styleUrls: ['./home.component.css']
})
export class HomeComponent {}
