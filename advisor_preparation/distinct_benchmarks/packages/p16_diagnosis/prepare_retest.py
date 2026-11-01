#!/usr/bin/env python3
"""Create and verify an additive source release, or audit actual-source subset results."""
import argparse,hashlib,json,shutil,tomllib
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def hashes(p):return {str(f.relative_to(p)):sha(f) for f in p.rglob('*') if f.is_file()}
def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--proof',type=Path,required=True);p.add_argument('--isolated',type=Path);p.add_argument('--audit',action='store_true');p.add_argument('--attempt',type=Path,required=True);a=p.parse_args()
 proof=json.loads(a.proof.read_text());assert proof['status']=='PASS' and proof['formula_permutation_checks']==4096 and proof['maximum_corrected_quspin_error']<=1e-10
 own=Path(__file__).parent;receipt=json.loads((own/'production_patch_receipt.json').read_text());old=a.root/'source/renyi_corrected';new=a.root/'source/renyi_disjoint_corrected'
 if a.audit:
  prod=tomllib.loads((a.attempt/'production_subsets.toml').read_text());assert prod['source_sha256']==receipt['candidate_sha256']
  expected={tuple(r['sites']):r['quspin_purity'] for r in proof['all_subset_purities']};rows=[dict(r,quspin_purity=expected[tuple(r['sites'])],absolute_error=abs(r['purity']-expected[tuple(r['sites'])])) for r in prod['rows']];assert len(rows)==16
  maximum=max(r['absolute_error'] for r in rows);report={'status':'PASS' if maximum<=1e-8 else 'FAIL','maximum_absolute_error':maximum,'rows':rows,'actual_source_sha256':prod['source_sha256'],'production_result_sha256':sha(a.attempt/'production_subsets.toml'),'proof_sha256':sha(a.proof),'scope':'All16 spatial subsets from actual repaired user source vs independent QuSpin; diagnostic support for P16, not extra ideas'}
  with (a.attempt/'production_subsets_comparison.json').open('x') as f:json.dump(report,f,indent=2)
  print(json.dumps({'status':report['status'],'maximum_absolute_error':maximum}));assert report['status']=='PASS';return
 frozen={'user_ed_distinct.jl':'d84fb923d3b29f2b978ae79378e0a59192b856c16c8218856b04fa63a1c50178','compare.py':'6a4b6bcc42a798cff152d4a5372935dbf6afcfc7e43eb30b787d1d2396d1ce50','manifest.toml':'6d3fd97ee4e83be6914416d2b645c3cb237f614638c4dfbb509516fbc14e71fe'}
 for name,digest in frozen.items():assert sha(own.parent/name)==digest,name
 assert sha(a.root/'results/distinct_tranche01_attempt01/oracle/oracle.toml')=='709a1a8303ecde7803283b7392e8fa7844ad4752c7c0ead7f3525cc5f1a0828b'
 assert a.isolated is not None
 isolated=tomllib.loads(a.isolated.read_text());assert isolated['status']=='PASS' and isolated['source_sha256']==receipt['original_sha256'] and isolated['maximum_absolute_error']<=1e-8
 assert sha(old/'analysis/compute_s2_ed.jl')==receipt['original_sha256'];assert sha(own/'compute_s2_ed_v2.jl')==receipt['candidate_sha256'];assert sha(own/'compute_s2_ed_v2.patch')==receipt['patch_sha256']
 assert not new.exists(),'New source release already exists; preserve it'
 shutil.copytree(old,new,symlinks=True);shutil.copy2(own/'compute_s2_ed_v2.jl',new/'analysis/compute_s2_ed.jl')
 before,after=hashes(old),hashes(new);assert before.keys()==after.keys();changed=[k for k in before if before[k]!=after[k]];assert changed==['analysis/compute_s2_ed.jl'],changed
 assert sha(new/'analysis/compute_s2_ed.jl')==receipt['candidate_sha256']
 record={'schema':'disjoint_source_release_v1','old_root':str(old),'new_root':str(new),'changed_paths':changed,'old_hashes':before,'new_hashes':after,'proof_sha256':sha(a.proof),'isolated_result_sha256':sha(a.isolated),'patch_receipt':receipt}
 a.attempt.mkdir(parents=True,exist_ok=False)
 with (a.attempt/'source_release_receipt.json').open('x') as f:json.dump(record,f,indent=2)
 print(json.dumps({'new_root':str(new),'changed_paths':changed,'new_source_sha256':receipt['candidate_sha256']}))
if __name__=='__main__':main()
