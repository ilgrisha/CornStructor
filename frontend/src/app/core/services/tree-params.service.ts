// File: frontend/src/app/core/services/tree-params.service.ts
// Version: v0.2.0
/**
 * TreeParamsService
 * - Holds Globals & Levels for construction tree generation.
 * - Loads defaults from /api/params/defaults (locked schema).
 * - Import/Export JSON (Globals and Levels separately).
 * - Add/Insert/Delete levels with index reflow and default value seeding.
 */
import { Injectable, signal, WritableSignal } from '@angular/core';
import { HttpClient } from '@angular/common/http';

export type JsonObject = { [k: string]: any };

export interface TreeDefaults { globals: JsonObject; levels: JsonObject[]; }

@Injectable({ providedIn: 'root' })
export class TreeParamsService {
  constructor(private http: HttpClient) {}

  readonly globals: WritableSignal<JsonObject> = signal({});
  readonly levels: WritableSignal<JsonObject[]> = signal([]);
  private defaults?: TreeDefaults;

  /** Fetch defaults from backend and populate working copies. */
  loadDefaults() {
    return this.http.get<TreeDefaults>('/api/params/defaults').subscribe(def => {
      this.defaults = def;
      this.globals.set(JSON.parse(JSON.stringify(def.globals)));
      this.levels.set(JSON.parse(JSON.stringify(def.levels)));
      this.reflowLevels();
    });
  }

  /** Export helpers. */
  exportGlobals(): string { return JSON.stringify(this.globals(), null, 2); }
  exportLevels(): string  { return JSON.stringify(this.levels(),  null, 2); }

  /** Import (applies immediately) with locked schema enforcement. */
  importGlobals(jsonText: string) {
    const data = JSON.parse(jsonText);
    if (!this.defaults) throw new Error('Defaults not loaded');
    this.globals.set(this.lockSchema(this.defaults.globals, data));
  }
  importLevels(jsonText: string) {
    const arr = JSON.parse(jsonText);
    if (!Array.isArray(arr)) throw new Error('Levels import must be an array');
    if (!this.defaults) throw new Error('Defaults not loaded');
    const tmpl = this.defaults.levels[0] ?? {};
    const locked = arr.map((it: any) => this.lockSchema(tmpl, it));
    this.levels.set(locked);
    this.reflowLevels();
  }

  /** Locked schema: only keys from template; preserve order. */
  private lockSchema(template: JsonObject, input: JsonObject): JsonObject {
    const out: JsonObject = {};
    for (const k of Object.keys(template)) {
      out[k] = (input && Object.prototype.hasOwnProperty.call(input, k)) ? input[k] : template[k];
    }
    return out;
  }

  /** Level operations */
  addLevel() {
    if (!this.defaults || this.defaults.levels.length === 0) return;
    const def = JSON.parse(JSON.stringify(this.defaults.levels[0]));
    const arr = [...this.levels()];
    arr.push(def);
    this.levels.set(arr);
    this.reflowLevels();
  }

  insertLevel(index: number) {
    if (!this.defaults || this.defaults.levels.length === 0) return;
    const def = JSON.parse(JSON.stringify(this.defaults.levels[0]));
    const arr = [...this.levels()];
    arr.splice(index, 0, def);
    this.levels.set(arr);
    this.reflowLevels();
  }

  deleteLevel(index: number) {
    const arr = [...this.levels()];
    arr.splice(index, 1);
    this.levels.set(arr);
    this.reflowLevels();
  }

  /** Ensure sequential 0..N for the `level` field if present. */
  reflowLevels() {
    const arr = [...this.levels()];
    arr.forEach((it, i) => {
      if (Object.prototype.hasOwnProperty.call(it, 'level')) it['level'] = i;
    });
    this.levels.set(arr);
  }
}
