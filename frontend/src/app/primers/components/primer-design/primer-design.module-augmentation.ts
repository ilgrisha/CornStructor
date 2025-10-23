// File: frontend/src/app/primers/components/primer-design/primer-design.module-augmentation.ts
// Version: v0.1.0
/**
 * Module augmentation: import SharedModule so the Log view is available.
 * Import this file in PrimersModule if your project structure separates SharedModule location.
 */

import { NgModule } from '@angular/core';
import { SharedModule } from '../../../shared/shared.module';

@NgModule({
  imports: [SharedModule],
  exports: [SharedModule]
})
export class PrimersSharedAugmentationModule {}
