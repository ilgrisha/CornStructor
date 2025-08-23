export interface FastaRecord { id: string; seq: string; desc?: string; }

export function parseFasta(text: string): FastaRecord[] {
  const lines = (text||'').split(/\r?\n/);
  const recs: FastaRecord[] = [];
  let id = ''; let desc = ''; let seq: string[] = [];
  const push = () => { if (id){ recs.push({ id, desc, seq: seq.join('').replace(/\s+/g,'').toUpperCase() }); } };
  for(const line of lines){
    if(line.startsWith('>')){
      push();
      const head = line.substring(1).trim();
      const sp = head.indexOf(' ');
      id = sp>0 ? head.substring(0, sp) : head;
      desc = sp>0 ? head.substring(sp+1) : '';
      seq = [];
    }else{
      seq.push(line.trim());
    }
  }
  push();
  return recs;
}
