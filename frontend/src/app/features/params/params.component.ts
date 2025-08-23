import { Component } from '@angular/core';
import { CardComponent } from '../../shared/ui/card/card.component';
import { AnalysisService } from '../../core/services/analysis.service';
import { NgIf } from '@angular/common';
import { FeatureKey } from '../../core/models/analysis';

@Component({
  selector: 'app-params',
  standalone: true,
  imports: [CardComponent, NgIf],
  templateUrl: './params.component.html',
  styleUrls: ['./params.component.css'],
})
export class ParamsComponent {
  constructor(public a: AnalysisService) {}
  p = this.a.params;
  tog = this.a.toggles;

  upd<K extends keyof ReturnType<typeof this.p>>(_key: K, _ev: Event) {
    // helper only for template type narrowing
  }
  updateNum<K extends keyof ReturnType<typeof this.p>>(key: K, ev: Event) {
    const val = Number((ev.target as HTMLInputElement).value);
    // ts-expect-error – generic numeric update
    this.a.updateParam(key, Number.isFinite(val) ? val : this.p()[key]);
  }
  toggle(key: FeatureKey, ev: Event) {
    const checked = (ev.target as HTMLInputElement).checked;
    this.a.updateToggle(key, checked);
  }
}
