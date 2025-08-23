/**
 * Path: frontend/src/app/features/home/home.component.ts
 * Version: v0.7.0
 *
 * Top-level screen layout:
 *  - Left: ParamsComponent (analysis parameters & toggles)
 *  - Right (top/bottom stack): SequenceInputComponent, ResultsComponent
 *  - Bottom (full width): RunComponent (Construction Tree Design)
 */

import { Component } from '@angular/core';
import { NgClass } from '@angular/common';
import { CardComponent } from '../../shared/ui/card/card.component';
import { SequenceInputComponent } from '../sequence-input/sequence-input.component';
import { ParamsComponent } from '../params/params.component';
import { ResultsComponent } from '../results/results.component';
import { RunComponent } from '../run/run.component';

@Component({
  selector: 'app-home',
  standalone: true,
  imports: [
    NgClass,
    CardComponent,
    SequenceInputComponent,
    ParamsComponent,
    ResultsComponent,
    RunComponent
  ],
  templateUrl: './home.component.html',
  styleUrls: ['./home.component.css'],
})
export class HomeComponent {}
