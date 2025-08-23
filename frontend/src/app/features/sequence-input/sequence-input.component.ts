/**
 * Path: frontend/src/app/features/sequence-input/sequence-input.component.ts
 * Version: v1.6.0
 *
 * Sequence editor card: textarea + FASTA upload + caret position.
 * Fixes:
 *  - Use AnalysisService public setters (setSequence, setCaret)
 *  - Read invalid bases from features().invalid via a computed
 */
import { Component, ElementRef, ViewChild, computed } from '@angular/core';
import { CardComponent } from '../../shared/ui/card/card.component';
import { NgIf } from '@angular/common';
import { AnalysisService } from '../../core/services/analysis.service';

@Component({
  selector: 'app-sequence-input',
  standalone: true,
  imports: [CardComponent, NgIf],
  templateUrl: './sequence-input.component.html',
  styleUrls: ['./sequence-input.component.css'],
})
export class SequenceInputComponent {
  constructor(public a: AnalysisService) {}

  @ViewChild('ta') ta!: ElementRef<HTMLTextAreaElement>;

  // invalid intervals from analysis service
  invalid = computed(() => this.a.features().invalid);

  get seq(): string { return this.a.sequence().toUpperCase(); }

  onInput() {
    const el = this.ta?.nativeElement;
    const raw = (el?.value ?? '').toUpperCase();
    // keep only letters and whitespace; analysis will mark non-ACGT in features
    const cleaned = raw.replace(/\s+/g, '');
    this.a.setSequence(cleaned);
    this.updateCaret();
  }

  onCaretEvent() { this.updateCaret(); }

  private updateCaret() {
    const pos = this.ta?.nativeElement?.selectionStart ?? null;
    this.a.setCaret(pos);
  }

  onFile(ev: Event) {
    const input = ev.target as HTMLInputElement;
    const file = input.files && input.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      const text = String(reader.result || '');
      // Parse first FASTA record if present; else treat as plain sequence
      let seq = '';
      if (text.trim().startsWith('>')) {
        const lines = text.split(/\r?\n/);
        for (const line of lines) {
          if (!line || line.startsWith('>')) continue;
          seq += line.trim();
        }
      } else {
        seq = text;
      }
      seq = seq.toUpperCase().replace(/\s+/g, '');
      this.a.setSequence(seq);
      // place caret at end
      queueMicrotask(() => {
        if (this.ta?.nativeElement) {
          const el = this.ta.nativeElement;
          el.value = seq;
          try {
            el.setSelectionRange(seq.length, seq.length);
          } catch {}
          this.a.setCaret(seq.length);
        }
      });
    };
    reader.readAsText(file);
    // reset input so selecting the same file re-triggers
    input.value = '';
  }
}
