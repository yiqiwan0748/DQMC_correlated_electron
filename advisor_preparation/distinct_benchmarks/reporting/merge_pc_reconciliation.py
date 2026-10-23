#!/usr/bin/env python3
"""Build the v12 scientific leaf from the immutable public39 base and pinned public receipts.

This is read-only analysis of existing results. It neither runs scientific solvers
nor submits jobs. Use the four packet replay.py entry points to reconstruct the
raw mathematics; their recorded audit decisions are independently gated here.
"""
from pathlib import Path,PurePosixPath
import argparse,collections,hashlib,json,tarfile,copy,sys
BASE_SHA='abee1b7af4263979880d24afc679d1553567e2141caf16c41888c7c46cf2e228'
CHANGED={'X09':'PASS','X13':'PASS','X16':'PASS','P21':'UNRESOLVED','X01':'UNRESOLVED'}
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def safe(path):
 p=PurePosixPath(path)
 if p.is_absolute() or '..' in p.parts or '\\' in path or ':' in path:raise ValueError('Unsafe path')
 return p

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo-root',type=Path,default=HERE.parents[2]);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 root=a.repo_root.resolve();db=root/'advisor_preparation/distinct_benchmarks';packet=db/'pc_continuation_20261023';basepath=db/'reporting/reports/ground_initial_v10/ledger_100.json'
 if sha(basepath)!=BASE_SHA:raise ValueError('Wrong public39 base')
 base=load(basepath);assert base['counts']=={'PASS':39,'FAIL':0,'UNRESOLVED':1,'NOT_RUN':60}
 manifest=load(packet/'packet_manifest.json');checked=[]
 for row in manifest['files']:
  q=root.joinpath(*safe(row['path']).parts)
  if not q.is_file() or q.is_symlink() or q.stat().st_size!=row['bytes'] or sha(q)!=row['sha256']:raise ValueError('Input closure differs: '+row['path'])
  checked.append(row)
 # Preserve all existing receipt files and check every prior leaf input/analysis binding.
 base_inputs=[]
 for row in base['input_manifest'].values():
  q=basepath.parent/row['evidence_copy'];assert q.is_file() and q.stat().st_size==row['bytes'] and sha(q)==row['sha256'];base_inputs.append(row['evidence_copy'])
 for path,digest in base['analysis_receipts'].items():assert sha(basepath.parent/path)==digest
 # Validate complete regular-text archive membership, not merely archive hash labels.
 archive_checks=[]
 for mf in sorted(packet.glob('*/*manifest.json')):
  m=load(mf)
  if m.get('schema')!='scientific_text_archive_v1':continue
  arc=mf.parent/m['archive'];assert sha(arc)==m['archive_sha256'] and arc.stat().st_size==m['archive_bytes']
  expected={r['path']:r for r in m['files']};assert len(expected)==m['file_count'] and len({k.casefold() for k in expected})==len(expected)
  seen=set()
  with tarfile.open(arc,'r:gz') as t:
   for member in t:
    safe(member.name);assert member.isfile() and member.name not in seen and member.name in expected
    data=t.extractfile(member).read();row=expected[member.name];assert len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256']
    data.decode('utf-8');assert b'\0' not in data;seen.add(member.name)
  assert seen==set(expected);archive_checks.append({'path':str(arc.relative_to(root)),'sha256':sha(arc),'regular_text_members':len(seen)})
 assert len(archive_checks)==6 and sum(r['regular_text_members'] for r in archive_checks)==2158
 # Integrity PASS is never substituted for scientific PASS.
 x09=load(packet/'X09/public_replay_receipt.json');assert x09['status']=='PASS' and x09['scientific_status']=={'X09':'PASS','X12':'PASS'} and x09['new_distinct_passes']==1 and x09['check_count']==56
 p21=load(packet/'P21/public_replay_receipt.json');assert p21['audit_status']=='PASS' and p21['scientific_status']=='UNRESOLVED' and p21['scientific_completions']==0 and p21['all_1408_real_and_704_complex_statistics_agree']
 ps=load(packet/'P21/audit_summary.json');fg=ps['failure_gates'];assert len(fg['real_ineligible'])==2 and fg['complex_precision_failures']==376 and fg['complex_rectangles_not_contained']==170 and fg['complex_rectangles_wholly_outside']==0
 xx=load(packet/'X13_X16/public_replay_receipt.json');assert xx['status']=='PASS' and xx['scientific_status']=={'X13':'PASS','X16':'PASS'} and xx['new_distinct_ideas']==2 and xx['pilot_scientific_completions']==0
 assert len(xx['rows'])==3 and all(r['status']=='PASS' and r['frozen_comparator_gates_and_labels_exact'] for r in xx['rows']) and [r['job'] for r in xx['rows']]==['x13_6462017','x16pilot_6462342','x16prod_6462471']
 x01=load(packet/'X01/public_replay_receipt.json');assert x01['audit_status']=='PASS' and x01['scientific_status']=='UNRESOLVED' and x01['new_scientific_completions']==0 and x01['selected_stage']==2 and x01['check_count']==4857 and x01['failure_count']==0
 xs=load(packet/'X01/recorded_details/summary.json');assert xs['stages']['2']['status']=='UNRESOLVED' and xs['stages']['3']['no_primary_reconstruction'] and xs['stages']['3']['no_ledger_promotion']
 bad=[v for v in xs['stages']['2']['endpoint_precision'] if not v['eligible']];assert len(bad)==1 and bad[0]['name']=='S2_quartic_continuum' and bad[0]['halfwidth']>bad[0]['maximum']==.005
 selected=load(packet/'selected_rows.json');updates={r['id']:r for r in selected['rows']};assert set(updates)==set(CHANGED)
 old={r['id']:r for r in base['rows']};assert len(old)==100 and len({r['id'] for r in base['rows']})==100
 defining=('id','category','title','scientific_question','catalogue_model','catalogue_ensemble','source_urls')
 for ident,r in updates.items():
  assert r['status']==CHANGED[ident] and all(r[k]==old[ident][k] for k in defining)
  assert r['previous_public_row_sha256']==hashlib.sha256(canon(old[ident])).hexdigest()
  assert not any(k in r for k in ('administrative_submission_history','verifier_metadata_history','current_attempt','excluded_pilot_attempt'))
 rows=[copy.deepcopy(updates.get(r['id'],r)) for r in base['rows']]
 preserved={r['id']:hashlib.sha256(canon(r)).hexdigest() for r in base['rows'] if r['id'] not in CHANGED}
 assert len(preserved)==95 and all(canon(r)==canon(old[r['id']]) for r in rows if r['id'] not in CHANGED)
 prior_pass=[r['id'] for r in base['rows'] if r['status']=='PASS'];assert len(prior_pass)==39 and all(canon(next(r for r in rows if r['id']==i))==canon(old[i]) for i in prior_pass)
 counts={k:sum(r['status']==k for r in rows) for k in ('PASS','FAIL','UNRESOLVED','NOT_RUN')};assert counts=={'PASS':42,'FAIL':0,'UNRESOLVED':2,'NOT_RUN':56}
 cats={cat:{k:sum(r['category']==cat and r['status']==k for r in rows) for k in counts} for cat in sorted({r['category'] for r in rows})}
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 report={'schema':'public_scientific_pc_reconciliation_leaf_v12','generated_utc':manifest['frozen_utc'],'tool_sha256':sha(__file__),'counts':counts,'category_counts':cats,'counting_rule':base['counting_rule'],'scope_limit':'Three new finite-fixture ideas pass: X09 precision certification, X13 charged visibility and X16 conditional binding. P21 and X01 retain inadequate precision/mixing. Pilots, repeated X12, setup-only PySCF and cancelled/incomplete attempts add no completions. Newly added measurement adapters are disclosed; no whole-package or thermodynamic certification is implied.','base_ledger_sha256':BASE_SHA,'base_ledger':'../ground_initial_v10/ledger_100.json','base_preservation':{'all95unaffected_row_objects_unchanged':True,'all39prior_PASS_rows_unchanged':True,'canonical_preserved_row_sha256':preserved},'changed_ids':sorted(CHANGED),'import_provenance_only':{k:selected[k] for k in ('original_PC_ledger_sha256','private_import_zip_sha256')},'public_input_manifest':'../../../pc_continuation_20261023/packet_manifest.json','packet_manifest_sha256':sha(packet/'packet_manifest.json'),'family_replays':{'X09':'../../../pc_continuation_20261023/X09/public_replay_receipt.json','P21':'../../../pc_continuation_20261023/P21/public_replay_receipt.json','X13_X16':'../../../pc_continuation_20261023/X13_X16/public_replay_receipt.json','X01':'../../../pc_continuation_20261023/X01/public_replay_receipt.json'},'rows':rows}
 # The relative paths above start from the report directory.
 for key in ['public_input_manifest']:
  assert (db/'pc_continuation_20261023/packet_manifest.json').is_file()
 for family in report['family_replays']:assert (packet/family/'public_replay_receipt.json').is_file()
 write(out/'ledger_100.json',report)
 write(out/'audit.json',{'schema':'public_scientific_leaf_closure_audit_v12','status':'PASS','ledger_sha256':sha(out/'ledger_100.json'),'counts':counts,'changed_ids':sorted(CHANGED),'preserved_row_objects':95,'preserved_prior_PASS_rows':39,'new_public_files_checked':len(checked),'prior_leaf_input_files_checked':len(base_inputs),'prior_leaf_analysis_receipts_checked':len(base['analysis_receipts']),'scientific_text_archives':archive_checks,'archive_regular_text_members':sum(r['regular_text_members'] for r in archive_checks),'scientific_runs_executed':False,'private_import_required':False,'statistical_gates_changed':False,'source_bytes_changed':False})
 (out/'README.md').write_text('''# Reconciled scientific checkpoint:23 October2026\n\n**42 PASS,2 UNRESOLVED and56 NOT_RUN** across exactly100 distinct catalogue IDs. The [ledger](ledger_100.json) retains all95 unaffected row objects and all39 prior passing rows from the immutable [39-question checkpoint](../ground_initial_v10/README.md). Only X09,X13,X16,P21 andX01 are updated.\n\nX09 now passes the original finite four-site ground-certificate gates with the frozen sixteen-sweep precision release; all earlier failed gates remain in the prior leaf. X13 passes complete charged thresholds/operator visibility over128 channels. X16 passes eight finite-fixture binding intervals under its explicitly conditional binary64 numerical budget. Its pilot remains excluded; no thermodynamic pairing claim follows.\n\nP21 remains UNRESOLVED: two real components fail the individual-chain ESS gate,376 complex intervals miss precision and170 rectangles are not wholly contained in the equivalence disk. Its actual production protocol explicitly supersedes the older unexecuted Welch/pilot candidate. X01 remains UNRESOLVED at stage2 because quartic entropy halfwidth0.0051318999469 exceeds0.005; its incomplete final stage3 contributes no primary estimate or credit. All original failed analyses, repairs, partial outputs and cancellations remain scientific history.\n\nThe [scientific packets](../../../pc_continuation_20261023/README.md) contain byte-exact selected source/raw text in safe archives, complete member manifests, unchanged mathematical auditors, frozen comparators and fresh portable replay receipts. Original import/PC-ledger hashes are provenance only; no private ZIP, controller, CV or contact data is required. The [closure audit](audit.json) verifies complete archives and all prior receipt bindings. This leaf adds no new scientific execution.\n'''.replace('checkpoint:23','checkpoint: 23').replace('October2026','October 2026').replace('PASS,2','PASS, 2').replace('and56','and 56').replace('exactly100','exactly 100').replace('all95','all 95').replace('all39','all 39').replace('X09,X13,X16,P21 andX01','X09, X13, X16, P21 and X01').replace('over128','over 128').replace('gate,376','gate, 376').replace('and170','and 170').replace('stage2','stage 2').replace('halfwidth0.','halfwidth 0.').replace('exceeds0.','exceeds 0.').replace('stage3','stage 3'),encoding='utf-8')
 print(json.dumps({'status':'PASS','counts':counts,'ledger_sha256':sha(out/'ledger_100.json'),'public_files':len(checked),'archive_members':2158,'output':str(out)},indent=2))
if __name__=='__main__':main()
