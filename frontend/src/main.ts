/* ============================================================================
 * Path: frontend/src/main.ts
 * Version: v3.0.0
 * ==========================================================================*/
import { enableProdMode, importProvidersFrom } from '@angular/core';
import { bootstrapApplication } from '@angular/platform-browser';
import { provideHttpClient } from '@angular/common/http';
import { HomeComponent } from './app/features/home/home.component';

bootstrapApplication(HomeComponent, {
  providers: [
    provideHttpClient(),
  ]
}).catch(err => console.error(err));
