/* ============================================================================
 * Path: frontend/src/app/features/sequence-box/sequence-box.component.ts
 * Version: v2.5.0
 * ============================================================================
 * - Dynamic per-line width (fits container; perfect 1ch alignment)
 * - Single interactive box (type/paste A/C/G/T, FASTA upload)
 * - Visible blinking caret; keyboard nav & editing
 * - Mouse selection across multiple lines; copy/cut/paste/replace
 * - Smooth scroll to center on selected feature-range change
 * - Bars: faint “all ranges” (optional), solid selected-range
 * - Meta (length/selected/cursor/FASTA) shown under the box
 * ==========================================================================*/
import {
  Component, computed, effect, signal, WritableSignal,
  ViewChild, ElementRef, AfterViewInit, OnDestroy
} from '@angular/core';
import { CommonModule } from '@angular/common';
import { AnalysisService, FeatureRegion } from '../../core/services/analysis.service';

type Cell = { ch: string; cls: string; abs: number };
type Line = {
  start: number;
  cells: Cell[];
  ruler: string[];
  barAll: boolean[];
  barSel: boolean[];
};

@Component({
  selector: 'app-sequence-box',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './sequence-box.component.html',
  styleUrls: ['./sequence-box.component.css']
})
export class SequenceBoxComponent implements AfterViewInit, OnDestroy {
  constructor(public a: AnalysisService) {
    // re-run grid when deps change
    effect(() => { void this.lines(); });

    // center-scroll when a feature range is selected from the table
    effect(() => {
      const r = this.selRange();
      if (r) {
        // move caret to the beginning of the selected range
        this.cursorPos.set(r.start);
        this.clearSelection();
        this.scheduleScrollTo(r.start);
      }
    });
  }

  @ViewChild('viz', { static: true }) viz!: ElementRef<HTMLDivElement>;
  private ro?: ResizeObserver;
  private winHandler = () => this.measureAndSetLineWidth();

  // ---- layout & state
  lineWidth: WritableSignal<number> = signal(80);
  hideAllWhenSelected: WritableSignal<boolean> = signal(true);

  // caret + selection (absolute indices; selectionEnd is exclusive)
  cursorPos = signal<number>(0);
  private selecting = false;
  private selectionAnchor: number | null = null;
  private selectionFocus: number | null = null;

  // invalid stats
  hasInvalid = computed(() => this.a.invalid().some(v => v));
  invalidCount = computed(() => this.a.invalid().reduce((s, v) => s + (v ? 1 : 0), 0));

  // FASTA meta
  private _lastFastaName: WritableSignal<string | null> = signal(null);
  private _lastFastaId: WritableSignal<string | null> = signal(null);
  lastFastaName = computed(() => this._lastFastaName());
  lastFastaId = computed(() => this._lastFastaId());

  // selected feature ranges from service
  selRanges = this.a.selectedRegions;
  selRange = this.a.selectedRange;

  // grid lines
  lines = computed<Line[]>(() =>
    this.buildLines(
      this.a.sequence().toUpperCase(),
      this.selRanges(),
      this.selRange(),
      this.lineWidth()
    )
  );

  // --- UI handlers
  onHideAllToggle(ev: Event) {
    const checked = (ev.target as HTMLInputElement).checked;
    this.hideAllWhenSelected.set(checked);
  }
  focusViz() { this.viz.nativeElement.focus(); }

  onFastaFile(ev: Event) {
    const file = (ev.target as HTMLInputElement).files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      const text = String(reader.result ?? '');
      // basic FASTA parse (first record)
      const m = /^>([^\n\r]*)[\r\n]+([\s\S]*)$/m.exec(text);
      let id: string | null = null;
      let seq = text;
      if (m) { id = m[1].trim(); seq = m[2]; }
      const cleaned = seq.replace(/[^ACGTacgt]/g, '').toUpperCase();
      this.a.sequence.set(cleaned);
      this._lastFastaName.set(file.name);
      this._lastFastaId.set(id);
      this.cursorPos.set(cleaned.length);
      this.clearSelection();
    };
    reader.readAsText(file);
    (ev.target as HTMLInputElement).value = '';
  }

  ngAfterViewInit(): void {
    this.measureAndSetLineWidth();
    this.ro = new ResizeObserver(() => this.measureAndSetLineWidth());
    this.ro.observe(this.viz.nativeElement);
    window.addEventListener('resize', this.winHandler);

    // selection end on mouseup anywhere
    window.addEventListener('mouseup', this.onGlobalMouseUp, true);

    // focus to enable typing immediately
    setTimeout(() => this.viz.nativeElement.focus());
  }
  ngOnDestroy(): void {
    this.ro?.disconnect();
    window.removeEventListener('resize', this.winHandler);
    window.removeEventListener('mouseup', this.onGlobalMouseUp, true);
  }

  // ----- selection helpers
  private getSelRange(): [number, number] | null {
    if (this.selectionAnchor == null || this.selectionFocus == null) return null;
    const a = this.selectionAnchor;
    const b = this.selectionFocus;
    return a <= b ? [a, b] : [b, a];
  }
  private clearSelection() {
    this.selectionAnchor = null;
    this.selectionFocus = null;
  }
  private isSelectedAbs(abs: number): boolean {
    const r = this.getSelRange();
    if (!r) return false;
    return abs >= r[0] && abs < r[1];
  }

  // mouse selection across lines
  onCellMouseDown(abs: number, e: MouseEvent) {
    e.preventDefault();
    this.viz.nativeElement.focus();
    if (e.shiftKey && this.selectionAnchor != null) {
      // extend from existing anchor
      this.selectionFocus = abs + 1;
    } else {
      this.selectionAnchor = abs;
      this.selectionFocus = abs + 1;
    }
    this.selecting = true;
    this.cursorPos.set(abs);
  }
  onCellMouseEnter(abs: number, e: MouseEvent) {
    if (!this.selecting) return;
    this.selectionFocus = abs + 1; // inclusive of hovered cell
  }
  private onGlobalMouseUp = () => {
    if (!this.selecting) return;
    this.selecting = false;
    const r = this.getSelRange();
    if (!r || r[0] === r[1]) this.clearSelection();
  };

  // caret keyboard & edit operations
  onVizKeydown(e: KeyboardEvent) {
    const key = e.key;
    const meta = e.metaKey || e.ctrlKey;

    // copy/cut
    if (meta && (key.toLowerCase() === 'c' || key.toLowerCase() === 'x')) {
      const sel = this.getSelRange();
      if (sel) {
        const s = this.a.sequence();
        const frag = s.slice(sel[0], sel[1]);
        navigator.clipboard?.writeText(frag).catch(() => {});
        if (key.toLowerCase() === 'x') {
          const ns = s.slice(0, sel[0]) + s.slice(sel[1]);
          this.a.sequence.set(ns);
          this.cursorPos.set(sel[0]);
          this.clearSelection();
        }
        e.preventDefault();
        return;
      }
    }

    // navigation & editing
    let s = this.a.sequence();
    let pos = this.cursorPos();
    const cols = this.lineWidth();

    const sel = this.getSelRange();

    const replaceSelection = (txt: string) => {
      if (!sel) return;
      const cleaned = txt.toUpperCase().replace(/[^ACGT]/g, '');
      const ns = s.slice(0, sel[0]) + cleaned + s.slice(sel[1]);
      this.a.sequence.set(ns);
      const newPos = sel[0] + cleaned.length;
      this.cursorPos.set(newPos);
      this.clearSelection();
      s = ns; pos = newPos;
    };

    const insertAtPos = (txt: string) => {
      const cleaned = txt.toUpperCase().replace(/[^ACGT]/g, '');
      if (!cleaned) return;
      const ns = s.slice(0, pos) + cleaned + s.slice(pos);
      this.a.sequence.set(ns);
      const newPos = pos + cleaned.length;
      this.cursorPos.set(newPos);
      s = ns; pos = newPos;
    };

    const backspace = () => {
      if (sel) {
        replaceSelection('');
        return;
      }
      if (pos <= 0) return;
      const ns = s.slice(0, pos - 1) + s.slice(pos);
      this.a.sequence.set(ns);
      this.cursorPos.set(pos - 1);
    };
    const del = () => {
      if (sel) {
        replaceSelection('');
        return;
      }
      if (pos >= s.length) return;
      const ns = s.slice(0, pos) + s.slice(pos + 1);
      this.a.sequence.set(ns);
    };

    // typing
    if (/^[acgtACGT]$/.test(key)) {
      if (sel) replaceSelection(key);
      else insertAtPos(key);
      e.preventDefault();
      return;
    }

    if (key === 'Backspace') { backspace(); e.preventDefault(); return; }
    if (key === 'Delete')    { del(); e.preventDefault(); return; }

    const moveCaret = (delta: number, keepSel = false) => {
      const np = Math.max(0, Math.min(s.length, pos + delta));
      if (keepSel) {
        if (this.selectionAnchor == null) this.selectionAnchor = pos;
        this.selectionFocus = np;
      } else {
        this.clearSelection();
      }
      this.cursorPos.set(np);
    };

    // arrow navigation (+ Shift to extend selection)
    if (key === 'ArrowLeft')  { moveCaret(-1, e.shiftKey); e.preventDefault(); return; }
    if (key === 'ArrowRight') { moveCaret(+1, e.shiftKey); e.preventDefault(); return; }
    if (key === 'ArrowUp')    { moveCaret(-cols, e.shiftKey); e.preventDefault(); return; }
    if (key === 'ArrowDown')  { moveCaret(+cols, e.shiftKey); e.preventDefault(); return; }
    if (key === 'Home')       { moveCaret(-(pos % cols), e.shiftKey); e.preventDefault(); return; }
    if (key === 'End')        { moveCaret(cols - (pos % cols), e.shiftKey); e.preventDefault(); return; }
    if (key === 'Enter')      { e.preventDefault(); return; }
  }

  // paste (replace selection if present)
  onVizPaste(e: ClipboardEvent) {
    const txt = e.clipboardData?.getData('text') ?? '';
    const cleaned = txt.replace(/[^ACGTacgt]/g, '').toUpperCase();
    if (!cleaned) { e.preventDefault(); return; }

    const sel = this.getSelRange();
    const s = this.a.sequence();
    const pos = this.cursorPos();

    if (sel) {
      const ns = s.slice(0, sel[0]) + cleaned + s.slice(sel[1]);
      this.a.sequence.set(ns);
      this.cursorPos.set(sel[0] + cleaned.length);
      this.clearSelection();
    } else {
      const ns = s.slice(0, pos) + cleaned + s.slice(pos);
      this.a.sequence.set(ns);
      this.cursorPos.set(pos + cleaned.length);
    }

    e.preventDefault();
  }

  // ----- scrolling to center a line
  private pendingScroll: number | null = null;
  private rafHandle: number | null = null;

  private scheduleScrollTo(absIndex: number) {
    this.pendingScroll = absIndex;
    if (this.rafHandle) cancelAnimationFrame(this.rafHandle);
    this.rafHandle = requestAnimationFrame(() => this.tryScrollPending());
  }

  private tryScrollPending() {
    if (this.pendingScroll == null) return;
    const viz = this.viz?.nativeElement;
    if (!viz) return;

    const cols = this.lineWidth();
    const rowIdx = Math.floor(this.pendingScroll / cols);

    const rows = Array.from(viz.querySelectorAll<HTMLElement>('.line'));
    if (!rows.length || !rows[rowIdx]) {
      // wait next frame if not yet rendered
      this.rafHandle = requestAnimationFrame(() => this.tryScrollPending());
      return;
    }
    const rowEl = rows[rowIdx];
    const rowTop = rowEl.offsetTop;
    const rowHeight = rowEl.offsetHeight;

    const targetTop = rowTop - (viz.clientHeight / 2 - rowHeight / 2);
    const clamped = Math.max(0, Math.min(viz.scrollHeight - viz.clientHeight, targetTop));
    viz.scrollTo({ top: clamped, behavior: 'smooth' });
    viz.focus();

    this.pendingScroll = null;
    if (this.rafHandle) { cancelAnimationFrame(this.rafHandle); this.rafHandle = null; }
  }

  // ----- layout computation
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

  // ----- building lines
  private cellClass(ch: string): string {
    if (ch === 'A') return 'nt nt-a';
    if (ch === 'C') return 'nt nt-c';
    if (ch === 'G') return 'nt nt-g';
    if (ch === 'T') return 'nt nt-t';
    return 'nt nt-x';
  }

  private buildRulerChars(start0: number, len: number): string[] {
    const arr = Array<string>(len).fill(' ');
    for (let col = 0; col < len; col++) {
      const abs1 = start0 + col + 1;
      if (abs1 % 10 === 0) {
        const s = String(abs1);
        const begin = col - (s.length - 1);
        for (let k = 0; k < s.length; k++) {
          const idx = begin + k;
          if (idx >= 0 && idx < len) arr[idx] = s[k];
        }
      }
    }
    return arr;
  }

  private paintBar(start0: number, len: number, ranges: FeatureRegion[] | null): boolean[] {
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

  private buildLines(
    seq: string, ranges: FeatureRegion[], selected: FeatureRegion | null, lineW: number
  ): Line[] {
    const out: Line[] = [];
    for (let i = 0; i < seq.length; i += lineW) {
      const chunk = seq.slice(i, i + lineW);
      const cells: Cell[] = Array.from(chunk, (ch, idx) => ({
        ch, cls: this.cellClass(ch), abs: i + idx
      }));
      out.push({
        start: i,
        cells,
        ruler: this.buildRulerChars(i, cells.length),
        barAll: this.paintBar(i, cells.length, ranges),
        barSel: this.paintSingle(i, cells.length, selected),
      });
    }
    return out;
  }

  // colors for bars (driven by selected kind)
  barColorAll = computed(() => {
    const k = this.a.selectedKind();
    return k ? this.a.featureColor(k) : 'transparent';
  });
  barColorSel = this.barColorAll;
}
