/* ============================================================================
 * Path: frontend/src/app/features/home/home.component.ts
 * Version: v2.3.0
 * ============================================================================
 * Layout update:
 * - Left (narrow column): Parameters + Detected features (Feature Picker)
 * - Right (wide column): Sequence (top) + Construction tree design (Run) below
 * ==========================================================================*/
import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';

import { ParamsComponent } from '../params/params.component';
import { FeaturePickerComponent } from '../feature-picker/feature-picker.component';
import { SequenceBoxComponent } from '../sequence-box/sequence-box.component';
import { RunComponent } from '../run/run.component';

@Component({
  selector: 'app-home',
  standalone: true,
  imports: [
    CommonModule,
    ParamsComponent,
    FeaturePickerComponent,
    SequenceBoxComponent,
    RunComponent,
  ],
  templateUrl: './home.component.html',
  styleUrls: ['./home.component.css'],
})
export class HomeComponent {}
