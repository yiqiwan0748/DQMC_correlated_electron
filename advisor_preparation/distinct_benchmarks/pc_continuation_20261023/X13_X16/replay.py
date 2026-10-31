#!/usr/bin/env python3
"""Portable unchanged bit-Fock/raw-math and frozen-comparator replays, no solver."""
import argparse,hashlib,importlib.util,json,math,os,subprocess,sys,tarfile
from pathlib import Path
from evidence_io import extract_verified,sha
HERE=Path(__file__).resolve().parent
PINS={'audit_x13_charged_spectral_v3.py':'b3b69b4e4fbcc1b009ef0b24ec9eb45303264ac6f792607a9c8772f2212bd80c','audit_x16_pair_binding_v3.py':'435b77332158a79172db1c6d7fc168ec3f290834ce915abb0d692e068baa731c'}
def module(name):
 p=HERE/'independent'/name
 if sha(p)!=PINS[name]:raise ValueError('Independent math source differs')
 spec=importlib.util.spec_from_file_location(name.replace('.','_'),p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def same(a,b):
 if isinstance(a,dict):
  if not isinstance(b,dict) or set(a)!=set(b):raise ValueError('Key grid differs')
  for k in a:
   if k=='input_sha256':
    if sorted(a[k].values())!=sorted(b[k].values()):raise ValueError('Input hash multiset differs')
   else:same(a[k],b[k])
 elif isinstance(a,list):
  if not isinstance(b,list) or len(a)!=len(b):raise ValueError('Component grid differs')
  for x,y in zip(a,b):same(x,y)
 elif isinstance(a,bool) or not isinstance(a,(float,int)):
  if a!=b:raise ValueError('Label/gate differs')
 elif not math.isfinite(a) or not math.isfinite(b) or abs(a-b)>64*sys.float_info.epsilon*max(1,abs(a),abs(b)):
  raise ValueError('Numeric replay beyond existing64epsilon allowance')
def plain(x):
 if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [plain(v) for v in x]
 if hasattr(x,'item'):return plain(x.item())
 return x
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
 manifest=json.loads((HERE/'data_manifest.json').read_text());selection=json.loads((HERE/'scientific_selection.json').read_text());assert selection['public_archive_members']==manifest['file_count']==113 and len(selection['excluded_files'])==33 and selection['scientific_source_and_raw_bytes_changed'] is False;assert not ({r['path'] for r in selection['excluded_files']} & {r['path'] for r in manifest['files']});raw=extract_verified(HERE/manifest['archive'],manifest,out/'raw');x13=module('audit_x13_charged_spectral_v3.py');x16=module('audit_x16_pair_binding_v3.py');rows=[]
 for job,key,m in [('x13_6462017','X13',x13),('x16pilot_6462342','X16',x16),('x16prod_6462471','X16',x16)]:
  root=raw/job;source=root/'source';freeze=m.load(source/'source_freeze.json');cfg=m.load(source/'manifest.toml');archive=HERE/'source_archives'/freeze['archive']
  if sha(archive)!=freeze['archive_sha256']:raise ValueError('Source archive differs')
  with tarfile.open(archive) as tf:
   members={}
   for item in tf:
    if not item.isfile() or item.name.startswith('/') or '..' in Path(item.name).parts:raise ValueError('Invalid source archive')
    name=Path(item.name).name
    if name in members:raise ValueError('Duplicate frozen source basename')
    members[name]=hashlib.sha256(tf.extractfile(item).read()).hexdigest()
  if members!=freeze['files'] or any(sha(source/n)!=h for n,h in members.items()):raise ValueError('Source closure differs')
  target=m.load(root/'results/target.toml');reference=m.load(root/'results/reference.json');execution=m.load(root/'results/execution_scientific.json');audit=m.Audit();freeze_sha=sha(source/'source_freeze.json')
  if key=='X13':
   t=m.audit_backend(audit,target,'native_ed_new_adapter',source,cfg,freeze_sha);r=m.audit_backend(audit,reference,'independent_quspin',source,cfg,freeze_sha)
   for i,(x,y) in enumerate(zip(target['site_annihilators'],reference['site_annihilators'])):audit.compare('cross canonical site annihilator '+str(i),m.maxabs(m.sparse(x)-m.sparse(y)),m.GATES['matrix_abs'])
   channels=m.compare_backends(audit,t,r);status=m.status(audit);findings={'cases':4,'complete_channels':len(channels),'hidden_lowest_channels':sum(c['hidden_lowest'] for c in channels)}
  else:
   t=m.audit_backend(audit,target,'target',cfg,freeze_sha,source,execution);r=m.audit_backend(audit,reference,'reference',cfg,freeze_sha,source,execution)
   bindings,grounds=m.cross_backends(audit,t,r);status,ready=m.scientific_status(audit,bindings);m.saved_comparison(audit,m.load(root/'results/comparison.json'),execution['attempt_mode'],status,ready,t,bindings,grounds);findings={'cases':4,'binding_intervals':bindings,'fixedN_ground_ensembles':len(grounds)}
  if status!='PASS' or not audit.integrity() or not all(c['passed'] for c in audit.numerical+audit.comparisons):raise ValueError('Independent numerical/scientific gates reject '+job)
  result={'schema':'portable_x13_x16_raw_math_replay_v1','job':job,'catalogue_id':key,'status':status,'counts':{'integrity_checks':len(audit.checks),'numerical_qualification_checks':len(audit.numerical),'cross_backend_checks':len(audit.comparisons)},'findings':findings,'integrity':audit.checks,'qualification':audit.numerical,'comparisons':audit.comparisons}
  (out/(job+'_math.json')).write_text(json.dumps(plain(result),indent=2,allow_nan=False)+'\n')
  fresh=out/(job+'_comparison.json');argv=[sys.executable,'-B',str(source/'compare.py'),'--manifest',str(source/'manifest.toml'),'--freeze-sha256',freeze_sha,'--target',str(root/'results/target.toml'),'--reference',str(root/'results/reference.json'),'--output',str(fresh)];env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
  if key=='X16':env['X16_ATTEMPT_MODE']=execution['attempt_mode']
  run=subprocess.run(argv,env=env,capture_output=True,text=True,timeout=180);(out/(job+'_comparator.stdout.txt')).write_text(run.stdout);(out/(job+'_comparator.stderr.txt')).write_text(run.stderr)
  if run.returncode:raise ValueError('Frozen comparator rejected '+job)
  same(m.load(root/'results/comparison.json'),m.load(fresh));rows.append({'job':job,'catalogue_id':key,'status':status,'pilot_excluded':execution.get('attempt_mode')=='pilot','counts':result['counts'],'frozen_comparator_gates_and_labels_exact':True})
 summary={'schema':'portable_x13_x16_scientific_replay_v1','status':'PASS','rows':rows,'scientific_status':{'X13':'PASS','X16':'PASS'},'new_distinct_ideas':2,'pilot_scientific_completions':0,'private_import_and_controller_checks_replayed':False,'excluded_controller_or_admission_bodies':33,'derived_scientific_execution_projections':3,'scientific_selection_manifest_sha256':sha(HERE/'scientific_selection.json'),'scientific_checks':'All original complete-grid independent bit-Fock, CAR, eigen-residual/projector, spectral and pair-density/interval gates rerun unchanged; explicit public path mapping only.','physics_executed':False}
 (out/'replay_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
