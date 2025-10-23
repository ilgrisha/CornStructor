// File: frontend/src/app/primers/primers.module.patch.ts
// Version: v0.1.0
/**
 * Patch: integrate SharedModule into PrimersModule without editing other imports above.
 * If you prefer, merge this into primers.module.ts directly.
 */

import { NgModule } from '@angular/core';
import { SharedModule } from '../shared/shared.module';

@NgModule({
  imports: [SharedModule],
  exports: [SharedModule]
})
export class PrimersModuleSharedPatch {}
