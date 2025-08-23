export interface StartDesignResponse { jobId: string; }
export interface ResultLink { kind: string; path: string; }
export interface DesignResultPayload {
  ok: boolean;
  done: boolean;
  jobId: string;
  outputDir: string;
  links: ResultLink[];
  fastaOut?: string;
}
