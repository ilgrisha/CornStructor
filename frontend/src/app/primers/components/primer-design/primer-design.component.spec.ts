// File: frontend/src/app/primers/components/primer-design/primer-design.component.spec.ts
// Version: v0.1.0
/**
 * Basic test: component creation and form defaults.
 */
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { PrimerDesignComponent } from './primer-design.component';
import { ReactiveFormsModule, FormsModule } from '@angular/forms';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { RouterTestingModule } from '@angular/router/testing';

describe('PrimerDesignComponent', () => {
  let component: PrimerDesignComponent;
  let fixture: ComponentFixture<PrimerDesignComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [PrimerDesignComponent],
      imports: [ReactiveFormsModule, FormsModule, HttpClientTestingModule, RouterTestingModule]
    }).compileComponents();

    fixture = TestBed.createComponent(PrimerDesignComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create and set default params', () => {
    expect(component).toBeTruthy();
    expect(component.form.value.primerLengthMin).toBe(18);
    expect(component.form.value.primerLengthMax).toBe(25);
  });
});
