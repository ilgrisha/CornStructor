/**
 * Path: frontend/src/app/features/sequence-view/sequence-view.component.ts
 * Version: v1.1.0
 *
 * Multi-line DNA visualization with nucleotide coloring (A/C/G/T)
 * and in-line highlight for the currently selected feature region.
 * Renders a “Selected region” panel under the sequence.
 */
import { Component, computed } from '@angular/core';
import { NgFor, NgIf } from '@angular/common';
import { AnalysisService, FeatureKind } from '../../core/services/analysis.service';

@Component({
  selector: 'app-sequence-view',
  standalone: true,
  imports: [NgFor, NgIf],
  templateUrl: './sequence-view.component.html',
  styleUrls: ['./sequence-view.component.css'],
})
export class SequenceViewComponent {
  constructor(public a: AnalysisService) {}

  // Wrap width per line
  wrap = 60;

  // Precompute line ranges
  lines = computed(() => {
    const s = this.a.sequence();
    const L = s.length;
    const out: { start: number; end: number }[] = [];
    for (let i = 0; i < L; i += this.wrap) {
      out.push({ start: i, end: Math.min(i + this.wrap, L) });
    }
    return out;
  });

  // Helpers
  coord(start: number, end: number) {
    return `${start + 1}–${end}`;
  }

  // Per-nt color classes
  ntClass(nt: string) {
    switch (nt) {
      case 'A': return 'ntA';
      case 'C': return 'ntC';
      case 'G': return 'ntG';
      case 'T': return 'ntT';
      default:  return 'ntX';
    }
  }

  // Is a position inside selected region?
  inSelected(pos: number): { on: boolean; color?: string } {
    const sel = this.a.currentRegion();
    if (!sel) return { on: false };
    if (pos >= sel.start && pos < sel.end) {
      const color = this.a.featureColors[sel.kind as FeatureKind] || '#ffd54f';
      return { on: true, color };
    }
    return { on: false };
  }

  // Build an array for a line to *ngFor each char with styles
  buildLine(start: number, end: number) {
    const seq = this.a.sequence();
    const row = [];
    for (let i = start; i < end; i++) {
      const ch = seq[i] ?? '';
      const sel = this.inSelected(i);
      row.push({
        ch,
        cls: this.ntClass(ch),
        style: sel.on ? {'background-color': this.tint(sel.color!, 0.25)} : null,
      });
    }
    return row;
  }

  // Make a soft translucent tint from a solid hex color
  tint(hex: string, alpha = 0.2) {
    // accept #RRGGBB
    const m = /^#?([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i.exec(hex);
    if (!m) return `rgba(255,213,79,${alpha})`;
    const r = parseInt(m[1], 16);
    const g = parseInt(m[2], 16);
    const b = parseInt(m[3], 16);
    return `rgba(${r}, ${g}, ${b}, ${alpha})`;
  }

  // Selected region snippet (under the sequence)
  selectedSnippet = computed(() => {
    const sel = this.a.currentRegion();
    if (!sel) return null;
    const seq = this.a.sequence();
    return {
      ...sel,
      text: seq.slice(sel.start, sel.end),
      len: sel.end - sel.start,
      color: this.a.featureColors[sel.kind],
    };
  });
}
