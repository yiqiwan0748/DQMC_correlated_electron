#!/usr/bin/env python3
"""Compare actual user-ED outputs to independent QuSpin references, per distinct ID."""
import argparse,hashlib,json,math,os,tomllib
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--oracle-dir',type=Path,required=True);ap.add_argument('--target',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();root=Path(__file__).parent;mfile=root/'manifest.toml';m=tomllib.loads(mfile.read_text());oracle=tomllib.loads((a.oracle_dir/'oracle.toml').read_text());oreceipt=json.loads((a.oracle_dir/'receipt.json').read_text());target=tomllib.loads(a.target.read_text());assert oreceipt['manifest_sha256']==target['metadata']['manifest_sha256']==sha(mfile)
 assert abs(oreceipt['ground_energy_uncentered']-target['metadata']['ground_energy_uncentered'])<=1e-8
 rows=[]
 for c in m['cases']:
  cid=c['id'];row={'id':cid,'scientific_question':c['observable'],'target':'user Hubbard ED state/Hamiltonian plus explicitly new analysis adapter','oracle':'independent QuSpin Hamiltonian/operators/partial trace'}
  if not c['ready']:row.update(status='NOT_RUN',reason='TeNPy/user observable adapter not implemented in this tranche')
  elif oracle.get(cid,{}).get('status')!='ORACLE_READY' or target.get(cid,{}).get('status')!='TARGET_READY':row.update(status='UNRESOLVED',oracle=oracle.get(cid),target=target.get(cid))
  else:
   x=target[cid]['values'];y=oracle[cid]['values'];row.update(target_values=x,oracle_values=y)
   if len(x)!=len(y):row.update(status='FAIL',reason=f'Observable/subspace dimension differs: target{len(x)},oracle{len(y)}')
   else:
    tolerances=[m['probability_tolerance'] if cid=='P01' and i<3 else m['absolute_tolerance'] for i in range(len(x))];errors=[abs(u-v) for u,v in zip(x,y)];good=all(math.isfinite(u) and math.isfinite(v) and e<=t for u,v,e,t in zip(x,y,errors,tolerances));row.update(status='PASS' if good else 'FAIL',component_count=len(x),maximum_absolute_error=max(errors,default=0),absolute_tolerances=tolerances,absolute_errors=errors)
  rows.append(row)
 result={'protocol_revision':m['protocol_revision'],'scope':'Seven distinct catalogue questions; each matrix/momentum/index/level arm is part of one question, not counted as an additional idea. Additional TeNPy IDs remain NOT_RUN. No DQMC validation claimed by this exact-state tranche.','manifest_sha256':sha(mfile),'oracle_receipt':oreceipt,'target_metadata':target['metadata'],'target_file_sha256':sha(a.target),'oracle_file_sha256':sha(a.oracle_dir/'oracle.toml'),'comparison_source_sha256':sha(__file__),'cases':rows,'counts':{s:sum(r['status']==s for r in rows) for s in ['PASS','FAIL','UNRESOLVED','NOT_RUN']}}
 assert not a.output.exists(),'Use a new result path';a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'counts':result['counts'],'rows':[(r['id'],r['status'],r.get('maximum_absolute_error')) for r in rows]}));raise SystemExit(1 if any(r['status']=='FAIL' for r in rows) else 0)
if __name__=='__main__':main()
