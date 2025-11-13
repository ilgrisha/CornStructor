// File: frontend/src/app/core/services/design.service.spec.ts
// Version: v0.2.0
/// <reference types="jasmine" />
/**
 * Minimal test to ensure start() normalizes jobId from both keys.
 */
import { TestBed } from '@angular/core/testing';
import { HttpClientTestingModule, HttpTestingController } from '@angular/common/http/testing';
import { DesignService } from './design.service';

describe('DesignService', () => {
  let svc: DesignService;
  let http: HttpTestingController;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule],
      providers: [DesignService],
    });
    svc = TestBed.inject(DesignService);
    http = TestBed.inject(HttpTestingController);
  });

  afterEach(() => http.verify());

  it('maps job_id -> jobId', (done: DoneFn) => {
    svc.start('ACGTACGTAC').subscribe(({ jobId }) => {
      expect(jobId).toBe('abc123');
      done();
    });
    const req = http.expectOne('/api/design/start');
    req.flush({ job_id: 'abc123' });
  });

  it('passes through jobId', (done: DoneFn) => {
    svc.start('ACGTACGTAC').subscribe(({ jobId }) => {
      expect(jobId).toBe('xyz789');
      done();
    });
    const req = http.expectOne('/api/design/start');
    req.flush({ jobId: 'xyz789' });
  });
});
