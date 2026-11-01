#!/usr/bin/env python3
"""Replay actual compressed X17 scientific data only; no producers or jobs."""
import argparse,hashlib,importlib.util,json,os,shutil,subprocess,sys,tarfile,tomllib,math
from pathlib import Path,PurePosixPath
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import numpy as np
from evidence_io import extract_verified
HERE=Path(__file__).resolve().parent
FREEZE='171c58736e7f0e59443b4429848827e35e81e0fa44c84ac2eb53dddd15e8fbcb'
ARCHIVE='5dfe6b1718a9af77fbfde134d09194c7594c1a144ce0d3bb6314fabaa3b13958'
AUDITORS={'pilot':'d3220adf8d731925e4a52ac00fb3175c37e6121c45d2ea15ade7e951e5f57c5c','production':'0dfa7af559a3e60c8c3167d8bf75ab07f34beda2f617fe2ece6623bd076ce4f1'}
ORIGINAL_AUDITS={'pilot':'59688ee4187434bea7f72e790015107bd95d47d593e3790b31410c9b6f0b76bd','production':'b617c4775c525c7e0328b6c9ed1c32f0e6212021af9daf885e0492092d348fb2'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):
 p=Path(p);return tomllib.loads(p.read_text()) if p.suffix=='.toml' else json.loads(p.read_text())
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def source(out):
 idx=load(HERE/'sources/index.json');freeze=load(HERE/'sources/source_freeze.json');arc=HERE/'sources'/idx['archive'];assert sha(HERE/'sources/source_freeze.json')==FREEZE and sha(arc)==ARCHIVE and arc.stat().st_size==idx['archive_bytes'];seen=set();out.mkdir()
 with tarfile.open(arc,'r:gz') as tar:
  for m in tar:
   assert m.isfile() and m.mode==0o644 and m.name.startswith(idx['prefix']);n=m.name[len(idx['prefix']):];p=PurePosixPath(n);assert not p.is_absolute() and '..' not in p.parts and '\\' not in n and ':' not in n and n not in seen and n in freeze['files'];data=tar.extractfile(m).read();data.decode('utf-8');assert b'\0' not in data and hashlib.sha256(data).hexdigest()==freeze['files'][n];q=out/n;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(data);seen.add(n)
 assert seen==set(freeze['files']);shutil.copyfile(HERE/'sources/source_freeze.json',out/'source_freeze.json');return out

def module(path):
 spec=importlib.util.spec_from_file_location('frozen_x17_pilot_auditor',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def pilot_fragment(raw,fresh):
 path=HERE/'independent/audit_pilot.py';assert sha(path)==AUDITORS['pilot'];text=path.read_text();definitions=text[text.index(' def need('):text.index(' try:')];fragment=text[text.index('  target=load('):text.index(' except Exception as ex:')]
 # One explicitly disclosed metadata lookup substitution. Original geometry,
 # fermion algebra, RDM, eigenstate, observable and resource math stays exact.
 original="before['runtime_critical_file_sha256']";assert fragment.count(original)==1
 rewritten=fragment.replace(original,"runtime_projection['original_runtime_critical_file_sha256']")
 code='def selected_original_pilot_math():\n checks=[]\n metrics=[]\n out={}\n source=p/"source"\n cfg=load(source/"manifest.toml")\n g=cfg["gates"]\n'+definitions+' if True:\n'+rewritten+' return dict(out=out,checks=checks,metrics=metrics)\n'
 mod=module(path);scope=dict(mod.__dict__);scope.update(p=raw,runtime_projection=load(HERE/'pilot/scientific_runtime_projection.json'));exec(compile(code,str(path)+'#public_metadata_selection_only','exec'),scope);result=scope['selected_original_pilot_math']();assert all(x['passed'] for x in result['checks']);result['fragment_original_sha256']=hashlib.sha256(fragment.encode()).hexdigest();result['fragment_public_sha256']=hashlib.sha256(rewritten.encode()).hexdigest();return result

def metric_differences(old,new):
 differences=[]
 def walk(x,y,path=()):
  assert type(x) is type(y),(path,type(x),type(y))
  if isinstance(x,dict):assert x.keys()==y.keys();[walk(v,y[k],path+(k,)) for k,v in x.items()]
  elif isinstance(x,list):assert len(x)==len(y);[walk(v,w,path+(i,)) for i,(v,w) in enumerate(zip(x,y))]
  elif isinstance(x,float):
   assert math.isfinite(x) and math.isfinite(y)
   if x!=y:differences.append(dict(path=list(path),original=x,replay=y,absolute_difference=abs(x-y)))
  else:assert x==y,path
 walk(old,new);return differences

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--role',choices=['pilot','production'],required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False);assert np.__version__=='2.4.6';role=a.role;ap=HERE/'independent'/('audit_'+role+'.py');assert sha(ap)==AUDITORS[role] and sha(HERE/role/'original_independent_audit.json')==ORIGINAL_AUDITS[role]
 raw=extract_verified(HERE/role/'data.tar.gz',load(HERE/role/'data_manifest.json'),out/'raw');src=source(raw/'source');execution=load(raw/'results/execution.json');assert execution['mode']==role and execution['status']=='PROCESSES_COMPLETE' and [x['exit_code'] for x in execution['processes']]==[0,0,0]
 frozen=raw/('frozen_pilot_audit' if role=='pilot' else 'frozen_full_replay');frozen.mkdir()
 if role=='pilot':cmd=[sys.executable,'-B',str(src/'audit_pilot.py'),'--pilot-dir',str(raw/'results'),'--output-dir',str(out/'new_pilot_analysis')]
 else:cmd=[sys.executable,'-B',str(src/'compare.py'),'--target',str(raw/'results/target.toml'),'--reference',str(raw/'results/reference.json'),'--output',str(frozen/'comparison.json')]
 r=subprocess.run(cmd,capture_output=True,text=True,timeout=300,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'));(out/'frozen_analysis.stdout.txt').write_text(r.stdout);(out/'frozen_analysis.stderr.txt').write_text(r.stderr);assert r.returncode==0,'Frozen arithmetic replay did not pass'
 if role=='pilot':
  for n in ['admission.json','comparison_replay.json']:shutil.copyfile(out/'new_pilot_analysis'/n,frozen/n)
 manifest=out/'public_scientific_input_manifest.json';files={str(q.relative_to(raw)):dict(sha256=sha(q),bytes=q.stat().st_size) for q in sorted(raw.rglob('*')) if q.is_file()};write(manifest,dict(schema='PUBLIC_X17_SCIENTIFIC_INPUT_MANIFEST_V1',files=files))
 if role=='production':
  cmd=[sys.executable,'-B',str(ap),'--receipt-root',str(raw),'--output',str(out/'public_independent_audit.json'),'--public-input-manifest',str(manifest),'--expected-public-manifest-sha256',sha(manifest)];r=subprocess.run(cmd,capture_output=True,text=True,timeout=300,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'));(out/'auditor.stdout.txt').write_text(r.stdout);(out/'auditor.stderr.txt').write_text(r.stderr);assert r.returncode==0
  audit=load(out/'public_independent_audit.json');assert audit['audit_status']==audit['integrity']=='PASS' and audit['scientific_completions']==1 and audit['catalogue_status']=='PASS'
 else:
  project=load(HERE/'pilot/scientific_runtime_projection.json');selection=load(HERE/'pilot/scientific_selection.json');assert project['derived_from_original_transfer_sha256']==selection['original_before_transfer_sha256'];result=pilot_fragment(raw,frozen);audit=dict(result['out'],schema='X17_PUBLIC_SELECTED_PILOT_MATH_V1',checks=result['checks'],check_count=len(result['checks']),case_metrics=result['metrics'],fragment_original_sha256=result['fragment_original_sha256'],fragment_public_sha256=result['fragment_public_sha256'],explicit_metadata_substitution="before['runtime_critical_file_sha256'] replaced by named derived scientific_runtime_projection.original_runtime_critical_file_sha256",remote_transport_checks_replayed=False);write(out/'public_independent_audit.json',audit)
 original=load(HERE/role/'original_independent_audit.json');diff=metric_differences(original['case_metrics'],audit['case_metrics']);saved=load(raw/'results/comparison.json');summary=dict(schema='X17_COMPRESSED_PUBLIC_REPLAY_V1',role=role,audit_status='PASS',catalogue_status=audit['catalogue_status'],scientific_completions=0 if role=='pilot' else 1,pilot_completion_credit=0,source_archive_sha256=ARCHIVE,source_freeze_sha256=FREEZE,original_audit_sha256=ORIGINAL_AUDITS[role],auditor_sha256=AUDITORS[role],data_archive_sha256=sha(HERE/role/'data.tar.gz'),original_execution_sha256=sha(raw/'results/execution.json'),public_independent_audit_sha256=sha(out/'public_independent_audit.json'),public_input_manifest_sha256=sha(manifest),independent_check_count=audit['check_count'],original_independent_check_count=original['check_count'],frozen_gate_count=sum(len(c['gates']) for c in saved['cases']),complete_case_count=len(saved['cases']),complete_spin_sector_count=5*len(saved['cases']),case_metrics=audit['case_metrics'],case_metric_differences=diff,case_metrics_bitwise_identical=not diff,all_frozen_gates_pass=all(g['passed'] for c in saved['cases'] for g in c['gates']),scientific_producer_rerun=False,postprocessing_eigendecomposition=True,remote_transport_checks_replayed=False)
 if role=='pilot':summary.update(projected_full_seconds=audit['projected_full_seconds'],projected_full_peak_rss_KiB=audit['projected_full_peak_rss_KiB'])
 write(out/'replay_summary.json',summary);print(json.dumps({k:summary[k] for k in ['role','audit_status','catalogue_status','scientific_completions','independent_check_count','frozen_gate_count','case_metrics_bitwise_identical']},indent=2));print('Metric differences',len(diff))
if __name__=='__main__':main()
