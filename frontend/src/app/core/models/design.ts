// File: frontend/src/app/core/models/design.ts
// Version: v0.2.0
export interface Design {
  id: number;
  created_at: string;
  updated_at: string;
  sequence: string;
  sequence_len: number;
  params_json?: string | null;
  tree_json?: string | null;
  name?: string | null;
}

export interface DesignByRun extends Design {
  job_id: string;
}
