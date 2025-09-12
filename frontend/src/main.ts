/* ============================================================================
 * Path: frontend/src/main.ts
 * Version: v3.1.0
 * ============================================================================
 * Bootstraps AppComponent so the root selector remains <app-root> in index.html.
 * AppComponent renders <app-home>, whose selector is 'app-home'.
 * ==========================================================================*/
import { bootstrapApplication } from '@angular/platform-browser';
import { provideHttpClient } from '@angular/common/http';

import { AppComponent } from './app/app.component';

bootstrapApplication(AppComponent, {
  providers: [provideHttpClient()],
}).catch(err => console.error(err));
