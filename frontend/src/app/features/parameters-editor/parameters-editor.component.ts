// File: frontend/src/app/features/parameters-editor/parameters-editor.component.ts
// Version: v0.5.1
/**
 * ParametersEditorComponent (standalone modal)
 *
 * v0.5.1
 * - Removed `typeof` usage from template (not supported pre-Angular 19).
 * - Added `isBoolean()` helper and switched boolean branch to use it.
 */
import { Component, Input, Output, EventEmitter, signal } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { TreeParamsService } from '../../core/services/tree-params.service';

type ErrMap = Record<string, string | null>;
type JsonObject = { [k: string]: any };

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
    return typeof val === 'number' || (typeof val === 'string' && /^-?\d+(\.\d+)?$/.test(val));
  }
  isBoolean(val: any): boolean {
    return typeof val === 'boolean';
  }
  isObjectLike(val: any): boolean {
    return val !== null && typeof val === 'object' && !Array.isArray(val);
  }
  toPrettyJson(val: any): string {
    try { return JSON.stringify(val, null, 2); } catch { return String(val); }
  }
  /** Basic clamping for percentages/ratios; otherwise just numeric parsing. */
  clampNumber(key: string, raw: any): any {
    if (raw === '' || raw === null || raw === undefined) return raw;
    const num = typeof raw === 'number' ? raw : parseFloat(raw);
    if (Number.isNaN(num)) return raw;
    if (/pct|percent|rate|ratio/i.test(key)) return Math.max(0, Math.min(100, num));
    return num;
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
      this.tp.importLevels(String(reader.result));
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
  exportGlobals() { this.download('globals.json', this.tp.exportGlobals()); }
  exportLevels()  { this.download('levels.json',  this.tp.exportLevels()); }

  // ---------- Updates: Globals ----------
  onGlobalChange(key: string, value: any) {
    const v = this.clampNumber(key, value);
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
    const v = this.clampNumber(field, value);
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
  onLevelFieldChange(idx: number, key: string, value: any) {
    const v = this.clampNumber(key, value);
    const arr = [...this.tp.levels()];
    const row = { ...arr[idx], [key]: v };
    arr[idx] = row;
    this.tp.levels.set(arr);
  }

  onLevelJsonChange(idx: number, key: string, text: string) {
    if (!this.levelJsonErrs[idx]) this.levelJsonErrs[idx] = {};
    try {
      const parsed = JSON.parse(text);
      this.levelJsonErrs[idx][key] = null;
      const arr = [...this.tp.levels()];
      const row = { ...arr[idx], [key]: parsed };
      arr[idx] = row;
      this.tp.levels.set(arr);
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
