// File: frontend/src/app/primers/components/parameters-editor/parameters-editor.component.ts
// Version: v0.1.0
/**
 * PrimersParametersEditorComponent
 * --------------------------------
 * Modal dialog for editing primer design parameters (mirrors Construction Tree editor UX).
 * - Loads current parameters via GET /api/v1/primers/parameters on open.
 * - Saves via PUT /api/v1/primers/parameters.
 */
import { Component, EventEmitter, Input, OnChanges, OnDestroy, OnInit, Output, SimpleChanges } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormBuilder, ReactiveFormsModule, Validators } from '@angular/forms';
import { PrimersService } from '../../services/primers.service';
import { PrimerDesignParameters } from '../../models/primer-params.model';
import { Subscription } from 'rxjs';

@Component({
  selector: 'app-primers-parameters-editor',
  standalone: true,
  imports: [CommonModule, ReactiveFormsModule],
  templateUrl: './parameters-editor.component.html',
  styleUrls: ['./parameters-editor.component.scss']
})
export class PrimersParametersEditorComponent implements OnInit, OnChanges, OnDestroy {
  @Input() open = false;
  @Output() close = new EventEmitter<void>();

  loading = false;
  saving = false;
  error: string | null = null;
  saved = false;

  private sub?: Subscription;

  form = this.fb.group({
    primerLengthMin: [18, [Validators.required, Validators.min(10), Validators.max(60)]],
    primerLengthMax: [25, [Validators.required, Validators.min(10), Validators.max(60)]],
    primerTmMin: [55, [Validators.required]],
    primerTmMax: [62, [Validators.required]],
    primerGCMin: [40, [Validators.required]],
    primerGCMax: [65, [Validators.required]],
    primerHomopolymerMax: [4, [Validators.required]],
    primerThreePrimeEndLength: [5, [Validators.required]],
    primerThreePrimePrimerMatchMax: [2, [Validators.required]],
    primerTargetMatchMax: [5, [Validators.required]],
    primerTargetMatchNumberMax: [2, [Validators.required]],
    primerSecondaryStructureDeltaGMin: [-9.0, [Validators.required]],
    primerTmDifferenceMax: [3.0, [Validators.required]],
    weights: this.fb.group({
      wTm: [1.0],
      wGC: [0.5],
      wDimer: [1.0],
      w3p: [1.25],
      wOff: [1.0],
    }),
    randomSeed: [null as number | null],
  });

  constructor(private fb: FormBuilder, private api: PrimersService) {}

  ngOnInit(): void {
    if (this.open) this.fetch();
  }

  ngOnChanges(changes: SimpleChanges): void {
    if (changes['open']?.currentValue === true) {
      this.fetch();
    }
  }

  ngOnDestroy(): void {
    this.sub?.unsubscribe();
  }

  private fetch(): void {
    this.loading = true;
    this.error = null;
    this.saved = false;
    this.sub?.unsubscribe();
    this.sub = this.api.getParameters().subscribe({
      next: (p) => {
        this.form.reset(p as any);
        this.loading = false;
      },
      error: (err) => {
        this.error = err?.error?.detail ?? 'Failed to load parameters.';
        this.loading = false;
      },
    });
  }

  save(): void {
    if (this.form.invalid) return;
    const payload = this.form.value as unknown as PrimerDesignParameters;
    this.saving = true;
    this.error = null;
    this.api.updateParameters(payload).subscribe({
      next: () => {
        this.saving = false;
        this.saved = true;
        // Keep modal open so the user sees “Saved”; they can close manually.
      },
      error: (err) => {
        this.error = err?.error?.detail ?? 'Failed to save parameters.';
        this.saving = false;
      },
    });
  }

  onClose(): void {
    this.close.emit();
  }
}
