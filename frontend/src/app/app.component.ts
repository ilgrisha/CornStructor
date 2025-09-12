// File: frontend/src/app/app.component.ts
// Version: v0.2.0
// Minor: Ensure root uses only <app-home>, with HomeComponent selector fixed.
import { Component } from '@angular/core';
import { HomeComponent } from './features/home/home.component';

@Component({
  selector: 'app-root',
  standalone: true,
  imports: [HomeComponent],
  template: `<div class="container"><app-home></app-home></div>`,
})
export class AppComponent {}
