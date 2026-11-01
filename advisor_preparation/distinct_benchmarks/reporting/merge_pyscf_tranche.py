#!/usr/bin/env python3
"""Conditionally add only fully audited P23/P24 to the immutable public42 leaf.

This tool verifies retained data and mathematical replay receipts only. It never
runs producers, submits jobs, rewrites prior leaves or credits pilots/repairs.
"""
import argparse,copy,hashlib,json,tarfile
from pathlib import Path,PurePosixPath
BASE_SHA='976272652f95647c133eed9c6261edeedfcf0833b38cf8bc44095398474d6338'
AUDIT_SHA='0ca893c030f775dd6d912f4a1b89e8a4e3cc3339e7180b23ccb0dba73e068ab6'
AUDITOR_SHA='4b33bd92bbef0a719f47a11d65a0db0a297b0374096efde97ccb62f43d54b0b7'
CHANGED={'P23','P24'}
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def safe(n):
 p=PurePosixPath(n)
 if not n or p.is_absolute() or '..' in p.parts or '\\' in n or ':' in n:raise ValueError('Unsafe member/input path')
 return p

def verify_manifest(root,m):
 checked=[]
 for r in m['files']:
  q=root.joinpath(*safe(r['path']).parts)
  if not q.is_file() or q.is_symlink() or q.stat().st_size!=r['bytes'] or sha(q)!=r['sha256']:raise ValueError('Input mismatch: '+r['path'])
  checked.append(r['path'])
 if len(checked)!=len(set(checked)):raise ValueError('Duplicate manifest path')
 return checked

def archive_members(archive,expected,prefix=''):
 seen=set();counts=0
 with tarfile.open(archive,'r:gz') as t:
  for member in t:
   safe(member.name);assert member.isfile() and member.mode==0o644
   assert member.name.startswith(prefix);name=member.name[len(prefix):];safe(name)
   assert name not in seen and name in expected
   data=t.extractfile(member).read();r=expected[name];assert hashlib.sha256(data).hexdigest()==r['sha256']
   if 'bytes' in r:assert len(data)==r['bytes']
   data.decode('utf-8');assert b'\0' not in data;seen.add(name);counts+=1
 assert seen==set(expected) and len({n.casefold() for n in seen})==len(seen)
 return counts

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo-root',type=Path,default=HERE.parents[2]);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.repo_root.resolve();db=root/'advisor_preparation/distinct_benchmarks';packet=db/'pyscf_external_20261023';basepath=db/'reporting/reports/pc_reconciliation_v12/ledger_100.json'
 assert sha(basepath)==BASE_SHA;base=load(basepath);assert base['counts']=={'PASS':42,'FAIL':0,'UNRESOLVED':2,'NOT_RUN':56}
 # Exact prior leaf and its input closure remain immutable and accessible.
 prior_manifest=(basepath.parent/base['public_input_manifest']).resolve();assert prior_manifest.is_relative_to(root) and sha(prior_manifest)==base['packet_manifest_sha256'];prior=verify_manifest(root,load(prior_manifest))
 original39=db/'reporting/reports/ground_initial_v10/ledger_100.json';older=load(original39);assert sha(original39)==base['base_ledger_sha256']
 prior_old=[]
 for r in older['input_manifest'].values():
  q=original39.parent/r['evidence_copy'];assert q.is_file() and sha(q)==r['sha256'] and q.stat().st_size==r['bytes'];prior_old.append(r['evidence_copy'])
 for n,h in older['analysis_receipts'].items():assert sha(original39.parent/n)==h
 manifest=load(packet/'packet_manifest.json');checked=verify_manifest(root,manifest);assert manifest['schema']=='PUBLIC_PYSCF_PACKET_V1'
 # Missing production or incomplete replay is a hard refusal, never pilot credit.
 replay=load(packet/'production/public_replay_receipt.json');audit=load(packet/'production/original_independent_audit.json');publicaudit=load(packet/'production/public_independent_audit.json');pilot=load(packet/'pilot_public_replay_receipt.json')
 assert sha(packet/'production/original_independent_audit.json')==AUDIT_SHA and sha(packet/'independent/audit_pyscf_production.py')==AUDITOR_SHA
 assert audit['audit_status']==audit['integrity']==publicaudit['audit_status']==publicaudit['integrity']==replay['audit_status']==replay['integrity']=='PASS'
 assert audit['catalogue_status']==publicaudit['catalogue_status']==replay['catalogue_status']=={'P23':'PASS','P24':'PASS'}
 assert all(x['complete_case_count']==7 and x['complete_manifold_count']==56 and x['pilot_completion_credit']==0 for x in [audit,publicaudit,replay])
 assert replay['all_frozen_gates_pass'] and replay['all_case_metrics_exact'] and not replay['scientific_producer_executed']
 assert audit['case_metrics']==publicaudit['case_metrics']==replay['case_metrics'] and replay['auditor_sha256']==AUDITOR_SHA and replay['original_audit_sha256']==AUDIT_SHA
 assert sha(packet/'production/public_independent_audit.json')==replay['fresh_public_independent_audit_sha256']
 assert sha(packet/'production/public_frozen_comparator_replay.json')==replay['fresh_comparator_sha256']
 assert sha(packet/'production/public_scientific_input_manifest.json')==replay['public_manifest_sha256']
 assert pilot['audit_status']=='PASS' and pilot['catalogue_status']=={'P23':'NOT_RUN','P24':'NOT_RUN'} and pilot['scientific_completions']==0 and pilot['original_execution_status']=='UNRESOLVED' and pilot['original_process_exits']==[0,0,1]
 # Verify every original source and raw member, including the failed first pilot.
 archives=[]
 for s in load(packet/'sources/index.json')['releases']:
  f=load(packet/s['freeze']);arc=packet/s['archive'];assert sha(arc)==s['archive_sha256'] and sha(packet/s['freeze'])==s['freeze_sha256'] and arc.stat().st_size==s['archive_bytes'];expected={n:{'sha256':h} for n,h in f['files'].items()};n=archive_members(arc,expected,f.get('archive_prefix',s['release']+'/'));archives.append({'archive':s['archive'],'members':n,'sha256':sha(arc)})
 data={}
 for name in ['failed_pilot','repaired_pilot','production']:
  m=load(packet/name/'data_manifest.json');arc=packet/name/m['archive'];assert m['schema']=='scientific_text_archive_v1' and sha(arc)==m['archive_sha256'] and arc.stat().st_size==m['archive_bytes'];n=archive_members(arc,{r['path']:r for r in m['files']});assert n==m['file_count'];archives.append({'archive':name+'/'+m['archive'],'members':n,'sha256':sha(arc)})
  if name=='production':
   with tarfile.open(arc,'r:gz') as t:
    for f in ['results/comparison.json','results/execution.json']:data[f]=json.loads(t.extractfile(f).read())
 # Bind every replay role to the actual retained source/data bytes, including
 # the new local replay output; no private import or invented execution role.
 release=next(s for s in load(packet/'sources/index.json')['releases'] if s['release']=='pyscf_two_orbital_production_v2')
 freeze=load(packet/release['freeze']);rawmf=load(packet/'production/data_manifest.json')
 expected_roles={'source/'+n:h for n,h in freeze['files'].items()}
 expected_roles['source/source_freeze.json']=release['freeze_sha256']
 expected_roles.update({r['path']:r['sha256'] for r in rawmf['files']})
 expected_roles['independent_replay/comparison.json']=sha(packet/'production/public_frozen_comparator_replay.json')
 pm=load(packet/'production/public_scientific_input_manifest.json');assert pm['schema']=='PUBLIC_PYSCF_SCIENTIFIC_INPUT_MANIFEST_V1'
 expected_public={n:h for n,h in expected_roles.items() if n.startswith(('source/','results/','independent_replay/'))}
 assert {n:r['sha256'] for n,r in pm['files'].items()}==expected_public
 assert publicaudit['input_sha256']==expected_roles
 assert sha(packet/'production/data.tar.gz')==replay['data_archive_sha256']
 for n in ['results/target.toml','results/reference.json','results/comparison.json','results/execution.json']:
  assert audit['input_sha256'][n]==expected_roles[n]
 comp=data['results/comparison.json'];execution=data['results/execution.json']
 assert comp['mode']=='production' and not comp['pilot_excluded'] and comp['exception'] is None and comp['completion_count']==2 and comp['numerically_eligible'] and all(c['passed'] for c in comp['checks']) and len(comp['manifold_receipts'])==56
 assert {r['id']:r['status'] for r in comp['rows']}=={'P23':'PASS','P24':'PASS'} and len(comp['checks'])==audit['frozen_comparator_checks']==publicaudit['frozen_comparator_checks']==replay['frozen_comparator_checks']
 assert execution['mode']=='production' and execution['status']=='PROCESSES_COMPLETE' and len(execution['processes'])==3 and [r['exit_code'] for r in execution['processes']]==[0,0,0] and execution['cpu_limit']==1
 for f,k in [('comparison.json','retained_comparison_sha256'),('target.toml','retained_target_sha256'),('reference.json','retained_reference_sha256')]:assert execution['outputs'][f]['sha256']==replay[k]
 assert audit['execution_original_sha256']==publicaudit['execution_original_sha256']==replay['execution_original_sha256']
 # Preserve catalogue definitions and every other row object, including both unresolved histories.
 old={r['id']:r for r in base['rows']};assert len(old)==100 and all(old[k]['status']=='NOT_RUN' for k in CHANGED)
 updates={}
 shared_scope={'geometry':'Six open 2x2 two-orbital plaquette fixtures plus one mandatory open 2x1 directed-transition witness','ensemble':'Fixed (Nup,Ndown)=(1,1) or (2,2); uniform complete ground manifold within that sector','interaction':'Only intraorbital U; U=0,2,6 for plaquette and U=1 direction witness; no Hund/interorbital/pair-hopping terms','hopping':'Original plaquette t1=-1,t2=-1.75,t3=-0.85,t4=-0.65,h=0','physical_energy':'Native centered H restored by U*N/2-U*M/4; no chemical-potential term','reference':'Actual PySCF 2.14.0 FCI Hamiltonian and spin contractions, transition RDM and creation/destruction action checks; dense complete eigensystems','new_adapters':'Dense spectral decomposition, spin operator and directed orbital-transition measurements around existing native kinetic/Hamiltonian APIs','attribution':'Existing native package and PySCF operators; no newly invented iterative eigensolver or native DQMC dynamics claim','coverage':'Seven complete fixtures and eight complete separated manifolds each; 56 manifold records count as two catalogue ideas','limits':'Finite-cluster, selected-interaction validation only; no general package, full-Hund, large-size or thermodynamic certification'}
 for ident in sorted(CHANGED):
  r=copy.deepcopy(old[ident]);own=[c for c in comp['checks'] if c['phase']==ident];assert own
  r.update(status='PASS',reported_status='PASS',reason=('Complete independent native/PySCF many-electron spectra, invariant manifold projectors and exact spin labels pass the frozen gates.' if ident=='P23' else 'Complete directed orbital-transition strengths and reverse-direction controls agree with independent PySCF RDM/action calculations, including the asymmetric analytic witness.'),fixtures=7,comparisons=len(own),metrics=copy.deepcopy(own),scope=copy.deepcopy(shared_scope),input_role='Frozen full production6465923; pilot6465868 and analysis repair6465910 excluded from credit',audit_issues=[],previous_public_row_sha256=hashlib.sha256(canon(old[ident])).hexdigest(),evidence={'packet_manifest':'../../../pyscf_external_20261023/packet_manifest.json','production_source_archive_sha256':audit['source_archive_sha256'],'production_source_freeze_sha256':audit['source_freeze_sha256'],'original_independent_audit_sha256':AUDIT_SHA,'independent_auditor_sha256':AUDITOR_SHA,'independent_check_count':audit['check_count'],'public_replay_check_count':publicaudit['check_count'],'full_fixture_count':7,'full_manifold_count':56},scientific_history=[{'job':'6465868','mode':'pilot','original_workflow':'UNRESOLVED','process_exits':[0,0,1],'event':'Both science producers succeeded; comparator JSON serialization failed. Original raw/source/traceback preserved.','credit':0},{'job':'6465910','mode':'analysis-only repair','event':'Strict NumPy Boolean serialization repair; no scientific producer rerun, no model/formula/gate change. Original failed comparator cost plus repair cost retained.','credit':0},{'job':'6465923','mode':'full production','event':'All seven fixtures, 56 complete low-energy manifolds, frozen full comparator and independent Fock-algebra/public-payload replay pass.','catalogue_credit':1}])
  updates[ident]=r
 rows=[updates.get(r['id'],copy.deepcopy(r)) for r in base['rows']];preserved={r['id']:hashlib.sha256(canon(r)).hexdigest() for r in base['rows'] if r['id'] not in CHANGED};assert len(preserved)==98
 assert all(canon(r)==canon(old[r['id']]) for r in rows if r['id'] not in CHANGED)
 assert all(canon(next(r for r in rows if r['id']==i))==canon(old[i]) for i in old if old[i]['status'] in ['PASS','UNRESOLVED'])
 counts={s:sum(r['status']==s for r in rows) for s in ['PASS','FAIL','UNRESOLVED','NOT_RUN']};assert counts=={'PASS':44,'FAIL':0,'UNRESOLVED':2,'NOT_RUN':54}
 cats={cat:{s:sum(r['category']==cat and r['status']==s for r in rows) for s in counts} for cat in sorted({r['category'] for r in rows})}
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 report=dict(schema='public_pyscf_finite_cluster_leaf_v13',generated_utc=manifest['frozen_utc'],tool_sha256=sha(__file__),counts=counts,category_counts=cats,counting_rule=base['counting_rule'],scope_limit='Only P23/P24 gain credit from complete finite-fixture production. Failed pilot and serialization repair are preserved with zero credit; prior X01/P21 remain unresolved unchanged. No full-package, native DQMC dynamics, full-Hund or thermodynamic certification.',base_ledger_sha256=BASE_SHA,base_ledger='../pc_reconciliation_v12/ledger_100.json',base_preservation=dict(all98unaffected_row_objects_unchanged=True,all42prior_PASS_rows_unchanged=True,both_prior_UNRESOLVED_rows_unchanged=True,canonical_preserved_row_sha256=preserved),changed_ids=sorted(CHANGED),public_input_manifest='../../../pyscf_external_20261023/packet_manifest.json',packet_manifest_sha256=sha(packet/'packet_manifest.json'),family_replays={'pilot_no_credit':'../../../pyscf_external_20261023/pilot_public_replay_receipt.json','full_production':'../../../pyscf_external_20261023/production/public_replay_receipt.json'},rows=rows)
 write(out/'ledger_100.json',report)
 write(out/'audit.json',dict(schema='public_pyscf_leaf_closure_audit_v13',status='PASS',ledger_sha256=sha(out/'ledger_100.json'),counts=counts,changed_ids=sorted(CHANGED),preserved_row_objects=98,preserved_prior_PASS_rows=42,preserved_prior_UNRESOLVED_rows=2,public_input_files_checked=len(checked),previous_leaf_input_files_checked=len(prior),original39_evidence_copies_checked=len(prior_old),original39_analysis_receipts_checked=len(older['analysis_receipts']),archives=archives,archive_regular_text_members=sum(r['members'] for r in archives),scientific_runs_executed=False,source_bytes_changed=False,scientific_gates_changed=False,pilot_repair_completion_credit=0))
 (out/'README.md').write_text('''# Independent PySCF finite-cluster checkpoint\n\n**44 PASS, 2 UNRESOLVED and 54 NOT_RUN** across exactly 100 distinct catalogue IDs. Only P23 and P24 change from the immutable [42-question checkpoint](../pc_reconciliation_v12/README.md); all 98 other row objects, all 42 previous passing rows and both unresolved histories are preserved exactly.\n\nSeven full-production fixtures compare the existing native two-orbital Hamiltonian with actual PySCF FCI contractions. P23 covers complete low-energy energy manifolds, invariant projectors and spin labels. P24 covers directed orbital-transition strengths, their reverse and sum rules, including a deliberately asymmetric analytic witness. The model includes intraorbital Hubbard U only. New dense eigensystem and observable adapters are explicitly attributed; the result concerns selected finite clusters, without a whole-package or thermodynamic claim.\n\nThe original pilot's comparator serialization failure remains in the [scientific packet](../../../pyscf_external_20261023/README.md), with both successful scientific producer outputs and its failed execution receipt. The analysis-only Boolean repair changes no scientific formula or gate and receives no credit. Full production and independent public-payload mathematical replay establish the two new statuses. The [closure audit](audit.json) verifies pinned source/raw archives and predecessor receipts.\n''')
 print(json.dumps({'status':'PASS','counts':counts,'ledger_sha256':sha(out/'ledger_100.json'),'public_files':len(checked),'preserved_rows':98,'archive_members':sum(r['members'] for r in archives)},indent=2))
if __name__=='__main__':main()
