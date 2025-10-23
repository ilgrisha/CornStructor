// File: frontend/src/main.ts
// Version: v3.2.0
/**
 * Standalone bootstrap + Router from AppRoutingModule.
 */
import { bootstrapApplication } from '@angular/platform-browser';
import { provideHttpClient } from '@angular/common/http';
import { importProvidersFrom } from '@angular/core';
import { AppComponent } from './app/app.component';
import { AppRoutingModule } from './app/app-routing.module';

bootstrapApplication(AppComponent, {
  providers: [
    provideHttpClient(),
    importProvidersFrom(AppRoutingModule),
  ],
}).catch(err => console.error(err));
