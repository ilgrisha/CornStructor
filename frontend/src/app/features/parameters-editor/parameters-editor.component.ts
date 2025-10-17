// File: frontend/src/app/features/parameters-editor/parameters-editor.component.ts
// Version: v0.6.0
/**
 * ParametersEditorComponent (standalone modal)
 *
 * v0.6.0
 * - Boolean toggle support in Levels and GA groups (not only TM).
 * - Alias support for levels keys:
 *     overlap_tm_min  <-> overlap_min_tm
 *     overlap_tm_max  <-> overlap_max_tm
 *   On import/edit/export, we mirror both aliases to minimize backend changes.
 * - Export now serializes levels with mirrored aliases to maximize compatibility.
 * - Safer number parsing and nullish handling.
 */
import { Component, Input, Output, EventEmitter, signal } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { TreeParamsService } from '../../core/services/tree-params.service';

type ErrMap = Record<string, string | null>;
type JsonObject = { [k: string]: any };
type JsonArray = any[];

@Component({
  selector: 'app-parameters-editor',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './parameters-editor.component.html',
  styleUrls: ['./parameters-editor.component.css']
})
export class ParametersEditorComponent {
  @Input() open = false;
  @Output() close = new EventEmitter<void>();

  /** Currently selected level index. */
  sel = signal(0);

  /** Expanded/collapsed state for the TM group box. */
  tmOpen = signal(true);

  /** Validation error states for JSON editors */
  globalJsonErrs: ErrMap = {};
  levelJsonErrs: Record<number, ErrMap> = {};

  constructor(public tp: TreeParamsService) {}

  // ---------- Type helpers ----------
  isPrimitive(val: any): boolean {
    const t = typeof val;
    return val === null || t === 'string' || t === 'number' || t === 'boolean';
  }
  isNumeric(val: any): boolean {
    // allow numbers and numeric strings
    return typeof val === 'number' || (typeof val === 'string' && /^-?\d+(\.\d+)?$/.test(val.trim()));
  }
  isBoolean(val: any): boolean {
    return typeof val === 'boolean';
  }
  isObjectLike(val: any): boolean {
    return val !== null && typeof val === 'object' && !Array.isArray(val);
  }
  isArray(val: any): boolean {
    return Array.isArray(val);
  }
  toPrettyJson(val: any): string {
    try { return JSON.stringify(val, null, 2); } catch { return String(val); }
  }
  /** Basic clamping for percentages/ratios; otherwise just numeric parsing. */
  clampNumber(key: string, raw: any): any {
    if (raw === '' || raw === null || raw === undefined) return raw;
    const num = typeof raw === 'number' ? raw : parseFloat(String(raw));
    if (Number.isNaN(num)) return raw;
    if (/pct|percent|rate|ratio/i.test(key)) return Math.max(0, Math.min(100, num));
    return num;
  }

  // ---------- Alias helpers (Levels) ----------
  /** Mirror known alias pairs into the same object without overwriting explicit values. */
  private mirrorLevelAliases(row: JsonObject): JsonObject {
    const out = { ...row };

    // tm min
    const vMinA = out['overlap_tm_min'];
    const vMinB = out['overlap_min_tm'];
    if (vMinA !== undefined && vMinB === undefined) out['overlap_min_tm'] = vMinA;
    if (vMinB !== undefined && vMinA === undefined) out['overlap_tm_min'] = vMinB;

    // tm max
    const vMaxA = out['overlap_tm_max'];
    const vMaxB = out['overlap_max_tm'];
    if (vMaxA !== undefined && vMaxB === undefined) out['overlap_max_tm'] = vMaxA;
    if (vMaxB !== undefined && vMaxA === undefined) out['overlap_tm_max'] = vMaxB;

    return out;
  }

  /** Normalize an entire levels array (apply alias mirroring to every row). */
  private normalizeLevelsArray(arr: JsonArray): JsonArray {
    return (arr || []).map((row: any) => this.mirrorLevelAliases(row ?? {}));
  }

  // ---------- Import / Export ----------
  onImportGlobals(ev: Event) {
    const file = (ev.target as HTMLInputElement).files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      this.tp.importGlobals(String(reader.result));
      this.globalJsonErrs = {};
      this.tmOpen.set(true);
    };
    reader.readAsText(file);
    (ev.target as HTMLInputElement).value = '';
  }

  onImportLevels(ev: Event) {
    const file = (ev.target as HTMLInputElement).files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      try {
        const raw = JSON.parse(String(reader.result));
        const arr = Array.isArray(raw) ? raw : [];
        const normalized = this.normalizeLevelsArray(arr);
        // push into TP as JSON string to reuse its import logic if needed, else set directly
        this.tp.importLevels(JSON.stringify(normalized));
      } catch {
        // if parsing fails, fall back to the existing importer (it may handle validation)
        this.tp.importLevels(String(reader.result));
      }
      this.levelJsonErrs = {};
      this.sel.set(0);
    };
    reader.readAsText(file);
    (ev.target as HTMLInputElement).value = '';
  }

  private download(filename: string, text: string) {
    const blob = new Blob([text], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
  }

  exportGlobals() {
    // pass-through; globals format already matches your new schema
    this.download('globals.json', this.tp.exportGlobals());
  }

  exportLevels()  {
    // Export with alias mirroring to satisfy both new and old keys.
    const levels = this.normalizeLevelsArray(this.tp.levels());
    this.download('levels.json', JSON.stringify(levels, null, 2));
  }

  // ---------- Updates: Globals ----------
  onGlobalChange(key: string, value: any) {
    // render booleans as toggles if you add them later; for now numbers/strings
    const v = this.isNumeric(value) ? this.clampNumber(key, value) : value;
    const curr = this.tp.globals();
    this.tp.globals.set({ ...curr, [key]: v });
  }

  onGlobalJsonChange(key: string, text: string) {
    try {
      const parsed = JSON.parse(text);
      this.globalJsonErrs[key] = null;
      const curr = this.tp.globals();
      this.tp.globals.set({ ...curr, [key]: parsed });
    } catch {
      this.globalJsonErrs[key] = 'Invalid JSON';
    }
  }

  /** Ensure a nested group exists (e.g., ga, tm). */
  ensureGroup(key: string) {
    const g = this.tp.globals();
    if (!g[key] || typeof g[key] !== 'object') {
      const updated = { ...g, [key]: {} as JsonObject };
      this.tp.globals.set(updated);
    }
  }

  /** Update a single nested field inside a globals group (e.g., ga.population). */
  onGlobalNestedChange(group: string, field: string, value: any) {
    this.ensureGroup(group);
    const g = { ...this.tp.globals() };
    const groupObj = { ...(g[group] as JsonObject) };
    const v = this.isNumeric(value) ? this.clampNumber(field, value) : value;
    groupObj[field] = v;
    g[group] = groupObj;
    this.tp.globals.set(g);
  }

  /** Boolean switch helper for nested fields. */
  onGlobalNestedToggle(group: string, field: string, checked: boolean) {
    this.ensureGroup(group);
    const g = { ...this.tp.globals() };
    const groupObj = { ...(g[group] as JsonObject) };
    groupObj[field] = checked;
    g[group] = groupObj;
    this.tp.globals.set(g);
  }

  // ---------- Updates: Levels ----------
  private setLevelRow(idx: number, row: JsonObject) {
    const arr = [...this.tp.levels()];
    arr[idx] = row;
    this.tp.levels.set(arr);
  }

  onLevelFieldChange(idx: number, key: string, value: any) {
    // Update one primitive field, then mirror aliases and save.
    const current = { ...(this.tp.levels()[idx] ?? {}) };

    let v: any = value;
    if (this.isBoolean(current[key])) {
      v = !!value; // trust the toggle
    } else if (this.isNumeric(value)) {
      v = this.clampNumber(key, value);
    }

    const updated = { ...current, [key]: v };
    const mirrored = this.mirrorLevelAliases(updated);
    this.setLevelRow(idx, mirrored);
  }

  onLevelJsonChange(idx: number, key: string, text: string) {
    if (!this.levelJsonErrs[idx]) this.levelJsonErrs[idx] = {};
    try {
      const parsed = JSON.parse(text);
      this.levelJsonErrs[idx][key] = null;
      const current = { ...(this.tp.levels()[idx] ?? {}) };
      const updated = { ...current, [key]: parsed };
      const mirrored = this.mirrorLevelAliases(updated);
      this.setLevelRow(idx, mirrored);
    } catch {
      this.levelJsonErrs[idx][key] = 'Invalid JSON';
    }
  }

  // ---------- Level ops ----------
  onAddLevel() {
    this.tp.addLevel();
    if (this.sel() >= this.tp.levels().length) this.sel.set(this.tp.levels().length - 1);
  }
  onInsertLevel(i: number, ev?: MouseEvent) {
    ev?.stopPropagation();
    this.tp.insertLevel(i);
    this.sel.set(i);
  }
  onDeleteLevel(i: number, ev?: MouseEvent) {
    ev?.stopPropagation();
    this.tp.deleteLevel(i);
    if (this.sel() >= this.tp.levels().length) this.sel.set(Math.max(0, this.tp.levels().length - 1));
  }
  onSelectLevel(i: number) { this.sel.set(i); }
  onReloadDefaults() {
    this.tp.loadDefaults();
    // Alias-mirror current levels after reload (in case defaults use either naming)
    const arr = this.normalizeLevelsArray(this.tp.levels());
    this.tp.levels.set(arr);
    this.globalJsonErrs = {};
    this.levelJsonErrs = {};
    this.sel.set(0);
    this.tmOpen.set(true);
  }

  // ---------- GA & TM helpers ----------
  gaKeys(): string[] {
    const g = this.tp.globals();
    return g?.ga && typeof g.ga === 'object' ? Object.keys(g.ga) : [];
  }
  tmKeys(): string[] {
    const g = this.tp.globals();
    return g?.tm && typeof g.tm === 'object' ? Object.keys(g.tm) : [];
  }
  gaValue(key: string): any { return (this.tp.globals().ga ?? {})[key]; }
  tmValue(key: string): any { return (this.tp.globals().tm ?? {})[key]; }
}
