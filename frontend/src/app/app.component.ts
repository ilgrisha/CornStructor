import { Component } from '@angular/core';
import { HomeComponent } from './features/home/home.component';

@Component({
  selector: 'app-root',
  standalone: true,
  imports: [HomeComponent],
  template: `<div class="container"><app-home /></div>`
})
export class AppComponent {}
