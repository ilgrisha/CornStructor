// File: frontend/src/app/app.component.ts
// Version: v0.3.0
/**
 * Root shell with a router outlet.
 * Keep it minimal to avoid interfering with feature layouts.
 */
import { Component } from '@angular/core';
import { RouterOutlet } from '@angular/router';

@Component({
  selector: 'app-root',
  standalone: true,
  imports: [RouterOutlet],
  template: `
    <div class="container">
      <router-outlet></router-outlet>
    </div>
  `,
})
export class AppComponent {}
