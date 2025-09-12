export const normalizeSeq = (raw: string): string =>
  (raw || '').toUpperCase().replace(/[^ACGTN]/g, (m) => (m === '\n' || m === '\r' || m === '\t' || m === ' ' ? '' : ''));

export const invalidPositions = (seq: string): number[] => {
  const bad: number[] = [];
  for (let i=0;i<seq.length;i++){
    const c = seq[i];
    if (!/[ACGTN]/.test(c)) bad.push(i);
  }
  return bad;
}

export const colorClassFor = (c: string) => `nt-${/[ACGT]/.test(c) ? c : 'N'}`;
