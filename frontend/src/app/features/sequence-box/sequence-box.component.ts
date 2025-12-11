// File: frontend/src/app/features/sequence-box/sequence-box.component.ts
// Version: v0.2.9

/* ============================================================================
 * Visible 1ch-aligned caret (blinking) rendered as an overlay grid with
 * len+1 columns so it can sit between bases and after the last base.
 * Custom multi-line selection/copy/paste, responsive per-line width,
 * ticks+ruler+bars perfect alignment.
 *
 * UPDATE (v0.2.9):
 * - After FASTA upload and sequence set, trigger stems recompute immediately:
 *     void this.a.refreshStems();
 *   (Your AnalysisService also auto-refreshes on edits via its debounced effect.)
 * ==========================================================================*/
import {
  Component, computed, effect, signal, WritableSignal,
  ViewChild, ElementRef, AfterViewInit, OnDestroy
} from '@angular/core';
import { CommonModule } from '@angular/common';
import { AnalysisService, FeatureRegion, DesignOverlayMode, DesignLevelStat, DesignFragmentRange, SequenceReferenceInfo } from '../../core/services/analysis.service';

type Cell = { ch: string; cls: string; abs: number };
type Line = {
  start: number;
  len: number;
  cells: Cell[];
  tick: boolean[];
  ruler: string[];
  barAll: boolean[];
  barSel: boolean[];
  caret: boolean[];   // len+1; true where caret should be drawn for this line
  senseSegments: OverlaySegment[];
  antisenseSegments: OverlaySegment[];
};

type OverlaySegment = {
  startCol: number;
  endCol: number;
  fragment: DesignFragmentRange;
};

@Component({
  selector: 'app-sequence-box',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './sequence-box.component.html',
  styleUrls: ['./sequence-box.component.css']
})
export class SequenceBoxComponent implements AfterViewInit, OnDestroy {
  constructor(public a: AnalysisService) { effect(() => { void this.lines(); }); }

  @ViewChild('viz', { static: true }) viz!: ElementRef<HTMLDivElement>;
  private ro?: ResizeObserver;
  private winHandler = () => this.measureAndSetLineWidth();

  lineWidth: WritableSignal<number> = signal(80);
  hideAllWhenSelected: WritableSignal<boolean> = signal(true);

  // caret / selection / validity
  cursorPos = signal<number>(0);
  hasInvalid = computed(() => this.a.invalid().some(v => v));
  invalidCount = computed(() => this.a.invalid().reduce((s, v) => s + (v ? 1 : 0), 0));

  private selStart: WritableSignal<number | null> = signal(null);
  private selEnd: WritableSignal<number | null> = signal(null);
  private anchor: WritableSignal<number | null> = signal(null);
  private dragging = false;

  selRangeAbs = computed<{ lo: number; hi: number } | null>(() => {
    const s = this.selStart(); const e = this.selEnd();
    if (s == null || e == null || s === e) return null;
    return s < e ? { lo: s, hi: e } : { lo: e, hi: s };
  });

  selRanges = this.a.selectedRegions;
  selRange = this.a.selectedRange;

  lines = computed<Line[]>(() => {
    const showDesign = this.a.designOverlayVisible();
    return this.buildLines(
      this.a.sequence().toUpperCase(),
      this.selRanges(),
      this.selRange(),
      this.lineWidth(),
      this.cursorPos(),   // include caret so we recompute when it moves
      showDesign ? this.a.designOverlayFragmentsSense() : [],
      showDesign ? this.a.designOverlayFragmentsAntisense() : []
    );
  });

  // ---------- selection helpers ----------
  isSelected(abs: number): boolean {
    const r = this.selRangeAbs();
    return !!r && abs >= r.lo && abs < r.hi;
  }
  clearSelection() {
    this.selStart.set(null); this.selEnd.set(null); this.anchor.set(null);
  }
  setCollapsedCaret(pos: number) {
    this.cursorPos.set(pos);
    this.clearSelection();
  }
  private setSelectionFromAnchor(toAbsInclusive: number) {
    const anc = this.anchor();
    if (anc == null) return;
    this.selStart.set(anc);
    this.selEnd.set(toAbsInclusive + 1); // half-open [lo, hi)
    const hi = Math.max(anc, toAbsInclusive) + 1;
    this.cursorPos.set(hi); // caret to end
  }

  // ---------- UI handlers ----------
  onHideAllToggle(ev: Event) {
    const checked = (ev.target as HTMLInputElement).checked;
    this.hideAllWhenSelected.set(checked);
  }

  overlayOptionLabel(stat: DesignLevelStat): string {
    const fragLabel = stat.fragments === 1 ? 'fragment' : 'fragments';
    const oligoLabel = stat.oligos === 1 ? 'oligo' : 'oligos';
    return `Level ${stat.level} — ${stat.fragments} ${fragLabel}, ${stat.oligos} ${oligoLabel}`;
  }

  onDesignLevelChange(ev: Event) {
    const value = (ev.target as HTMLSelectElement).value;
    if (!value) {
      this.a.setDesignOverlayLevel(null);
      return;
    }
    this.a.setDesignOverlayLevel(Number(value));
  }

  onDesignModeChange(ev: Event) {
    const value = (ev.target as HTMLSelectElement).value as DesignOverlayMode;
    this.a.setDesignOverlayMode(value);
  }

  clearDesignOverlaySelection() {
    this.a.setDesignOverlayLevel(null);
    this.a.selectDesignFragment(null);
  }

  onFragmentSegmentClick(fragment: DesignFragmentRange, ev: MouseEvent) {
    ev.stopPropagation();
    this.a.selectDesignFragment(fragment);
  }

  fragmentLabel(fragment: DesignFragmentRange): string {
    return fragment.isOligo ? 'Oligo' : 'Fragment';
  }

  formatTm(tm: number | null): string {
    if (tm == null || Number.isNaN(tm)) return '—';
    return `${tm.toFixed(1)} °C`;
  }

  segmentTooltip(fragment: DesignFragmentRange): string {
    return `${this.fragmentLabel(fragment)} ${fragment.id} (${fragment.length} bp)`;
  }

  referenceLabel(meta: SequenceReferenceInfo | null): string {
    if (!meta) return '';
    const prefix = meta.source === 'fasta' ? 'FASTA' : 'Reference';
    const baseLabel = meta.label?.trim().length
      ? meta.label.trim()
      : meta.source === 'run' && meta.id
        ? `Run ${meta.id}`
        : 'Loaded sequence';
    const idSuffix = meta.id ? ` (${meta.id})` : '';
    return `${prefix}: ${baseLabel}${idSuffix}`;
  }

  focusViz() { this.viz.nativeElement.focus(); }

  onFastaFile(ev: Event) {
    const file = (ev.target as HTMLInputElement).files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      const text = String(reader.result ?? '');
      const m = /^>([^\n\r]*)[\r\n]+([\s\S]*)$/m.exec(text);
      let id: string | null = null;
      let seq = text;
      if (m) { id = m[1].trim(); seq = m[2]; }
      const cleaned = seq.replace(/[^ACGTacgt]/g, '').toUpperCase();
      this.a.clearDesignOverlayData();
      this.a.sequence.set(cleaned);
      this.a.setSequenceReference({ source: 'fasta', label: file.name, id });
      this.cursorPos.set(cleaned.length);
      this.clearSelection();

      // NEW: trigger stems immediately after upload (debounced auto-refresh will also run)
      void this.a.refreshStems();
    };
    reader.readAsText(file);
    (ev.target as HTMLInputElement).value = '';
  }

  ngAfterViewInit(): void {
    this.measureAndSetLineWidth();
    this.ro = new ResizeObserver(() => this.measureAndSetLineWidth());
    this.ro.observe(this.viz.nativeElement);
    window.addEventListener('resize', this.winHandler);
    window.addEventListener('mouseup', this.onWindowMouseUp, true);
    setTimeout(() => this.viz.nativeElement.focus());
  }
  ngOnDestroy(): void {
    this.ro?.disconnect();
    window.removeEventListener('resize', this.winHandler);
    window.removeEventListener('mouseup', this.onWindowMouseUp, true);
  }

  private onWindowMouseUp = () => { this.dragging = false; };

  onCellMouseDown(abs: number, shift: boolean) {
    this.viz.nativeElement.focus();
    if (shift) {
      if (this.anchor() == null) this.anchor.set(this.cursorPos());
    } else {
      this.anchor.set(abs);
      this.selStart.set(abs);
      this.selEnd.set(abs);
    }
    this.dragging = true;
  }
  onCellMouseEnter(abs: number) {
    if (!this.dragging) return;
    this.setSelectionFromAnchor(abs);
  }
  onCellClick(abs: number) {
    if (!this.dragging) {
      this.setCollapsedCaret(abs);
      this.viz.nativeElement.focus();
    }
  }

  onVizKeydown(e: KeyboardEvent) {
    const key = e.key;
    const meta = e.metaKey || e.ctrlKey;
    if (meta && key.toLowerCase() === 'c') { this.copySelected(); e.preventDefault(); return; }

    const s = this.a.sequence();
    let pos = this.cursorPos();
    const cols = this.lineWidth();
    const sel = this.selRangeAbs();

    const replaceSel = (text: string) => {
      const cleaned = text.toUpperCase().replace(/[^ACGT]/g, '');
      const r = this.selRangeAbs();
      if (!r) return false;
      const ns = s.slice(0, r.lo) + cleaned + s.slice(r.hi);
      this.a.sequence.set(ns);
      const newPos = r.lo + cleaned.length;
      this.cursorPos.set(newPos);
      this.clearSelection();
      return true;
    };
    const insert = (text: string) => {
      const cleaned = text.toUpperCase().replace(/[^ACGT]/g, '');
      if (!cleaned) return;
      if (sel) { replaceSel(text); return; }
      const ns = s.slice(0, pos) + cleaned + s.slice(pos);
      this.a.sequence.set(ns);
      this.cursorPos.set(pos + cleaned.length);
    };
    const backspace = () => {
      if (sel) { replaceSel(''); return; }
      if (pos <= 0) return;
      const ns = s.slice(0, pos - 1) + s.slice(pos);
      this.a.sequence.set(ns);
      this.cursorPos.set(pos - 1);
    };
    const del = () => {
      if (sel) { replaceSel(''); return; }
      if (pos >= s.length) return;
      const ns = s.slice(0, pos) + s.slice(pos + 1);
      this.a.sequence.set(ns);
    };

    if (/^[acgtACGT]$/.test(key)) { insert(key); e.preventDefault(); return; }
    if (key === 'Backspace') { backspace(); e.preventDefault(); return; }
    if (key === 'Delete')    { del(); e.preventDefault(); return; }
    if (key === 'ArrowLeft') { this.clearSelection(); if (pos > 0) this.cursorPos.set(pos - 1); e.preventDefault(); return; }
    if (key === 'ArrowRight'){ this.clearSelection(); if (pos < s.length) this.cursorPos.set(pos + 1); e.preventDefault(); return; }
    if (key === 'ArrowUp')   { this.clearSelection(); this.cursorPos.set(Math.max(0, pos - cols)); e.preventDefault(); return; }
    if (key === 'ArrowDown') { this.clearSelection(); this.cursorPos.set(Math.min(s.length, pos + cols)); e.preventDefault(); return; }
    if (key === 'Home')      { this.clearSelection(); this.cursorPos.set(pos - (pos % cols)); e.preventDefault(); return; }
    if (key === 'End')       { this.clearSelection(); this.cursorPos.set(Math.min(s.length, pos - (pos % cols) + cols)); e.preventDefault(); return; }
    if (key === 'Enter')     { e.preventDefault(); return; }
  }

  onVizPaste(e: ClipboardEvent) {
    const txt = e.clipboardData?.getData('text') ?? '';
    const cleaned = txt.replace(/[^ACGTacgt]/g, '');
    if (cleaned) {
      const s = this.a.sequence();
      const r = this.selRangeAbs();
      if (r) {
        const ns = s.slice(0, r.lo) + cleaned.toUpperCase() + s.slice(r.hi);
        this.a.sequence.set(ns);
        this.cursorPos.set(r.lo + cleaned.length);
        this.clearSelection();
      } else {
        const pos = this.cursorPos();
        const ns = s.slice(0, pos) + cleaned.toUpperCase() + s.slice(pos);
        this.a.sequence.set(ns);
        this.cursorPos.set(pos + cleaned.length);
      }
      // No need to call refreshStems() here; AnalysisService debounced effect will run.
    }
    e.preventDefault();
  }

  onVizCopy(e: ClipboardEvent) {
    const r = this.selRangeAbs();
    if (!r) return;
    const text = this.a.sequence().slice(r.lo, r.hi);
    e.clipboardData?.setData('text/plain', text);
    e.preventDefault();
  }
  private async copySelected() {
    const r = this.selRangeAbs();
    if (!r) return;
    const text = this.a.sequence().slice(r.lo, r.hi);
    try { await navigator.clipboard.writeText(text); }
    catch {
      const ta = document.createElement('textarea');
      ta.value = text; ta.style.position = 'fixed'; ta.style.left = '-9999px';
      document.body.appendChild(ta); ta.select(); document.execCommand('copy'); document.body.removeChild(ta);
    }
  }

  // ---------- layout helpers ----------
  private measureAndSetLineWidth(): void {
    const host = this.viz?.nativeElement;
    if (!host) return;

    const cs = getComputedStyle(host);
    const padd = parseFloat(cs.paddingLeft) + parseFloat(cs.paddingRight);
    const border = parseFloat(cs.borderLeftWidth) + parseFloat(cs.borderRightWidth);
    const usable = Math.max(0, host.clientWidth - padd - border);

    const probe = document.createElement('span');
    probe.style.position = 'absolute';
    probe.style.visibility = 'hidden';
    probe.style.pointerEvents = 'none';
    probe.style.whiteSpace = 'pre';
    probe.style.fontFamily = 'ui-monospace, SFMono-Regular, Menlo, monospace';
    probe.style.fontSize = '13px';
    probe.textContent = '0'.repeat(200);
    document.body.appendChild(probe);
    const charW = probe.offsetWidth / 200;
    document.body.removeChild(probe);

    const cols = Math.max(20, Math.min(400, Math.floor(usable / charW)));
    if (cols !== this.lineWidth()) this.lineWidth.set(cols);
  }

  private cellClass(ch: string): string {
    if (ch === 'A') return 'nt nt-a';
    if (ch === 'C') return 'nt nt-c';
    if (ch === 'G') return 'nt nt-g';
    if (ch === 'T') return 'nt nt-t';
    return 'nt nt-x';
  }
  private buildTickFlags(start0: number, len: number): boolean[] {
    const arr = Array<boolean>(len).fill(false);
    for (let col = 0; col < len; col++) {
      const abs1 = start0 + col + 1;
      if (abs1 % 10 === 0) arr[col] = true;
    }
    return arr;
  }
  private buildRulerChars(start0: number, len: number): string[] {
    const arr = Array<string>(len).fill(' ');
    for (let col = 0; col < len; col++) {
      const abs1 = start0 + col + 1;
      if (abs1 % 10 === 0) {
        const s = String(abs1);
        const startCol = col - Math.floor((s.length - 1) / 2);
        for (let k = 0; k < s.length; k++) {
          const idx = startCol + k;
          if (idx >= 0 && idx < len) arr[idx] = s[k];
        }
      }
    }
    return arr;
  }
  private paintBar(
    start0: number,
    len: number,
    ranges: Array<{ start: number; end: number }> | null
  ): boolean[] {
    const out = new Array<boolean>(len).fill(false);
    if (!ranges || !ranges.length) return out;
    for (const r of ranges) {
      const s = Math.max(0, r.start - start0);
      const e = Math.min(len, r.end - start0);
      for (let i = s; i < e; i++) out[i] = true;
    }
    return out;
  }
  private paintSingle(start0: number, len: number, r: FeatureRegion | null): boolean[] {
    const out = new Array<boolean>(len).fill(false);
    if (!r) return out;
    const s = Math.max(0, r.start - start0);
    const e = Math.min(len, r.end - start0);
    for (let i = s; i < e; i++) out[i] = true;
    return out;
  }

  private buildOverlaySegments(
    start0: number,
    len: number,
    frags: DesignFragmentRange[]
  ): OverlaySegment[] {
    if (!frags.length) return [];
    const lineEnd = start0 + len;
    const segments: OverlaySegment[] = [];
    for (const frag of frags) {
      const segStart = Math.max(start0, frag.start);
      const segEnd = Math.min(lineEnd, frag.end);
      if (segStart >= segEnd) continue;
      segments.push({
        startCol: segStart - start0,
        endCol: segEnd - start0,
        fragment: frag,
      });
    }
    return segments;
  }

  private buildLines(
    seq: string,
    ranges: FeatureRegion[],
    selected: FeatureRegion | null,
    lineW: number,
    caretAbs: number,
    senseFragments: DesignFragmentRange[],
    antisenseFragments: DesignFragmentRange[]
  ): Line[] {
    const out: Line[] = [];
    for (let i = 0; i < seq.length; i += lineW) {
      const chunk = seq.slice(i, i + lineW);
      const cells: Cell[] = Array.from(chunk, (ch, idx) => ({ ch, cls: this.cellClass(ch), abs: i + idx }));
      const len = cells.length;

      const caretArr = new Array<boolean>(len + 1).fill(false);
      if (caretAbs >= i && caretAbs <= i + len) {
        caretArr[caretAbs - i] = true;
      }

      out.push({
        start: i,
        len,
        cells,
        tick: this.buildTickFlags(i, len),
        ruler: this.buildRulerChars(i, len),
        barAll: this.paintBar(i, len, ranges),
        barSel: this.paintSingle(i, len, selected),
        caret: caretArr,
        senseSegments: this.buildOverlaySegments(i, len, senseFragments),
        antisenseSegments: this.buildOverlaySegments(i, len, antisenseFragments),
      });
    }
    // Edge case: empty sequence → still no lines (caret handled by placeholder)
    return out;
  }

  // color for bars
  barColorAll = computed(() => this.a.selectedKind() ? this.a.featureColor(this.a.selectedKind()!) : 'transparent');
  barColorSel = this.barColorAll;
  readonly senseBarColor = '#38bdf8';
  readonly antisenseBarColor = '#ec4899';
}
