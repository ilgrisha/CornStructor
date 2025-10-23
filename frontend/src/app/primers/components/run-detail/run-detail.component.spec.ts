// File: frontend/src/app/primers/components/run-detail/run-detail.component.spec.ts
// Version: v0.1.0
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { RunDetailComponent } from './run-detail.component';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { RouterTestingModule } from '@angular/router/testing';

describe('RunDetailComponent', () => {
  let component: RunDetailComponent;
  let fixture: ComponentFixture<RunDetailComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [RunDetailComponent],
      imports: [HttpClientTestingModule, RouterTestingModule]
    }).compileComponents();

    fixture = TestBed.createComponent(RunDetailComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
