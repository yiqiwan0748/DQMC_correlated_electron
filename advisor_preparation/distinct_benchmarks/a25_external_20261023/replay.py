#!/usr/bin/env python3
"""Portable saved-data mathematical replay; no producer, pinned-runtime bypass or comparator rerun."""
import argparse,collections,hashlib,importlib.util,json,math,re,shutil,tarfile,tempfile,sys,platform
from pathlib import Path,PurePosixPath
HERE=Path(__file__).resolve().parent
AUDITS={'pilot':'7c23588be8ffe7dd7a470760805ba846615796c2f3040030d1c70a3f0298d95f','production':'a3396cb84f3553d38d90f5c2277ca01c19e6f6d4f68d600f0d9db12085404be4'}
AUDITORS={'pilot':'3cd20a52d5ba557b8a12c086addb5fa2bd4acd934f2ae70baf8facb6da08915d','production':'137e68f31f46d3fb0974357bc151bf1b2fc0a1501a40de580c4c3333d5210148'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def canon(v):return json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False)
def write(p,q):
 with Path(p).open('x') as f:json.dump(q,f,indent=2,allow_nan=False);f.write('\n')
def safe(n):
 p=PurePosixPath(n)
 if not n or p.is_absolute() or '..' in p.parts or '\\' in n or ':' in n:raise ValueError('Unsafe public member')
 return p

def extract(archive,dest,files,prefix=''):
 seen=set()
 with tarfile.open(archive,'r:gz') as tar:
  for m in tar:
   safe(m.name);assert m.isfile() and m.mode==0o644 and m.name.startswith(prefix)
   n=m.name[len(prefix):];safe(n);assert n in files and n not in seen;seen.add(n);b=tar.extractfile(m).read();spec=files[n];assert hashlib.sha256(b).hexdigest()==spec['sha256'] and ('bytes' not in spec or len(b)==spec['bytes']);b.decode('utf8');assert b'\0' not in b;p=dest.joinpath(*safe(n).parts);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
 assert seen==set(files)

def derived_receipt_equal(original,fresh,path=(),differences=None):
 if differences is None:differences=[]
 if type(original) is not type(fresh):raise ValueError('Derived receipt type changed '+repr(path))
 if isinstance(original,dict):
  if original.keys()!=fresh.keys():raise ValueError('Derived receipt keys changed '+repr(path))
  for k in original:derived_receipt_equal(original[k],fresh[k],path+(k,),differences)
 elif isinstance(original,list):
  if len(original)!=len(fresh):raise ValueError('Derived receipt length changed '+repr(path))
  for i,(x,y) in enumerate(zip(original,fresh)):derived_receipt_equal(x,y,path+(i,),differences)
 elif isinstance(original,float):
  allowance=64*sys.float_info.epsilon*max(1,abs(original),abs(fresh))
  if not math.isfinite(original) or not math.isfinite(fresh) or abs(original-fresh)>allowance:raise ValueError('Derived receipt arithmetic changed beyond declared portable allowance '+repr(path))
  if original!=fresh:differences.append(dict(path=list(path),original=original,fresh=fresh,absolute_difference=abs(original-fresh),allowed_portable_receipt_roundoff=allowance))
 elif original!=fresh:raise ValueError('Derived receipt nonfloat identity changed '+repr(path))
 return differences

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--role',choices=['pilot','production'],required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();role=a.role;d=HERE/role;out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);checks=[]
 def need(ok,name):
  checks.append(dict(name=name,passed=bool(ok)))
  if not ok:raise ValueError(name)
 sourceindex=read(HERE/'sources/index.json');specs={v['role']:v for v in sourceindex['releases']};orig=read(d/'original_independent_audit.json');need(sha(d/'original_independent_audit.json')==AUDITS[role] and orig['audit_status']==orig['integrity_status']=='PASS' and orig['auditor_sha256']==AUDITORS[role],'Actual original analysis report and unchanged original auditor identity');need(orig['check_count']==len(orig['checks']) and all(x['passed'] for x in orig['checks']),'All original recorded audit checks pass')
 data=read(d/'data_manifest.json');need(sha(d/'data.tar.gz')==data['archive_sha256'] and len(data['files'])==data['file_count'],'Actual compressed scientific data closure')
 with tempfile.TemporaryDirectory(prefix='a25-public-saved-math-') as temp:
  work=Path(temp);source=work/'source';source.mkdir();analysis=work/'analysis_source';analysis.mkdir()
  for name,dest in [('scientific',source),('pilot_analysis' if role=='pilot' else 'full_analysis',analysis)]:
   sp=specs[name];fr=read(HERE/sp['freeze']);need(sha(HERE/sp['archive'])==sp['archive_sha256']==fr['archive_sha256'] and sha(HERE/sp['freeze'])==sp['freeze_sha256'],'Exact frozen '+name+' release');extract(HERE/sp['archive'],dest,{n:{'sha256':h} for n,h in fr['files'].items()},sp['prefix']);shutil.copyfile(HERE/sp['freeze'],dest/'source_freeze.json')
  extract(d/'data.tar.gz',work,data['files']);raw=work/'results';returned=read(work/'analysis/audit.json');need(returned==orig,'Exact returned compute audit body retained in payload');old=read(raw/'comparison.json');need(old==read(d/'original_comparison.json') and old==read(work/'analysis/frozen_replay/comparison.json') and sha(d/'original_pinned_frozen_comparator_replay.json')==sha(work/'analysis/frozen_replay/comparison.json'),'Retained actual pinned compute comparison replay exact to original')
  replay=orig.get('frozen_full_comparator_replay') if role=='production' else next(x for x in orig['frozen_replays'] if x['label']=='comparison');need(replay['exit_code']==0 and not replay['timed_out'],'Original actual pinned compute replay succeeded')
  if role=='production':need(replay['comparison_sha256']==sha(d/'original_pinned_frozen_comparator_replay.json') and replay['exact_JSON_equal'],'Full retained replay bytes bound by original audit')
  execution=read(raw/'execution.json');need(execution['schema']=='A25_SCIENTIFIC_EXECUTION_V1' and execution['mode']==('pilot' if role=='pilot' else 'production') and execution['status']=='PROCESSES_COMPLETE' and execution['threads']==1 and execution['scientific_completions']==0 and all(x['exit_code']==0 and not x['timed_out'] for x in execution['processes']),'Complete original scientific execution, all role exits and zero runner credit')
  for n,v in execution['outputs'].items():need(sha(raw/n)==v['sha256'] and (raw/n).stat().st_size==v['bytes'],'Original execution-bound scientific output '+n)
  for n in ['native.toml','reference.json','comparison.json','execution.json']:
   bindings=[h for path,h in orig['input_sha256'].items() if path.endswith('/results/results/'+n)];need(len(bindings)==1 and bindings[0]==sha(raw/n),'Actual saved audit binds retained original '+n)
  modulefile=analysis/('audit_a25_pilot.py' if role=='pilot' else 'audit_a25_full.py');need(sha(modulefile)==AUDITORS[role],'Byte-unchanged original mathematical auditor module');module_spec=importlib.util.spec_from_file_location('unchanged_a25_saved_mathematics',modulefile);module=importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(module)
  plan=read(analysis/('audit_plan.json' if role=='pilot' else 'plan.json'));need(sha(analysis/('audit_plan.json' if role=='pilot' else 'plan.json'))==orig['audit_plan_sha256'],'Byte-unchanged original analysis plan');runtime=read(source/'runtime_pins.json');mapping={n:v['sha256'] for n,v in runtime['files'].items()};need(len(mapping)==119 and execution['runtime_before_sha256']==execution['runtime_after_sha256']==mapping,'Retained original runtime metadata maps; not a local runtime qualification')
  mathaudit=module.Audit();details=module.audit_math(source,raw,plan,mathaudit,mapping);differences=[]
  if role=='production':points,responses=details;derived_receipt_equal(orig['point_details'],points,('point_details',),differences);derived_receipt_equal(orig['response_details'],responses,('response_details',),differences);need(True,'Every independent derived point/response receipt obeys declared portable float comparison');need(len(points)==196 and len(responses)==56 and old['sector_solves']==980 and len(old['checks'])==6072,'All registered196points980sectors56responses6072gates complete')
  else:points=details;responses=[];derived_receipt_equal(orig['point_details'],points,('point_details',),differences);need(len(points)==14 and old['sector_solves']==98 and len(old['checks'])==518,'Pilot all14point metrics98sector timings518gates unchanged and excluded')
  need(bool(mathaudit.checks) and all(v['passed'] for v in mathaudit.checks),'All unchanged original mathematical gates pass')
  original_checks=collections.Counter(canon({k:v for k,v in row.items() if k!='value'}) for row in orig['checks']);new_checks=collections.Counter(canon({k:v for k,v in row.items() if k!='value'}) for row in mathaudit.checks);need(all(original_checks[k]>=v for k,v in new_checks.items()),'Every scientific check name/schema/type/decision/literal bound unchanged; fresh error observations satisfy original helpers')
  if role=='pilot':
   assessment=read(d/'original_pinned_resource_assessment.json');need(assessment==orig['resource_assessment']==read(work/'analysis/frozen_replay/assessment.json') and assessment['scientific_completions']==0 and assessment['catalogue_status']=='NOT_RUN','Retained actual pinned resource assessment is exact and zero-credit');assessment_replay=next(x for x in orig['frozen_replays'] if x['label']=='assessment');need(assessment_replay['exit_code']==0 and not assessment_replay['timed_out'],'Original pinned assessment replay succeeded');native=module.read(raw/'native.toml');warm=[s['elapsed_seconds'] for point in native['points'][1:] for s in point['spin_sectors']];procs={v['label']:v for v in execution['processes']};native_est=2*(procs['native']['wall_seconds']+980*max(warm));total=native_est+28*(procs['reference']['wall_seconds']+procs['comparison']['wall_seconds'])+30;memory=2*(max(v['max_rss_kib'] for v in procs.values())/1024+64);need(len(warm)==91 and assessment['native_forecast_seconds']==native_est and assessment['full_forecast_seconds']==total and assessment['full_forecast_mib']==memory,'Predeclared retained full91warm/cold resource formula reproduced exactly')
  else:assessment=None
  result=dict(schema='A25_PUBLIC_UNCHANGED_SAVED_MATH_REPLAY_V1',status='PASS',catalogue_status='NOT_RUN' if role=='pilot' else 'PASS',scientific_completions=0 if role=='pilot' else 1,pilot_completion_credit=0,role=role,original_independent_audit_sha256=AUDITS[role],original_unchanged_auditor_sha256=AUDITORS[role],original_check_count=orig['check_count'],public_scientific_check_count=len(mathaudit.checks),public_closure_check_count=len(checks),checks=checks,scientific_checks=mathaudit.checks,point_details=points,response_details=responses,retained_original_metric_bytes_exact=True,fresh_derived_point_and_response_metrics_exact=not differences,portable_receipt_roundoff_differences=differences,portable_receipt_policy_sha256=sha(HERE/'PORTABLE_RECEIPT_POLICY.json'),portable_python_version=platform.python_version(),original_pinned_compute_python_version='3.11.11',source_archive_sha256=specs['scientific']['archive_sha256'],source_freeze_sha256=specs['scientific']['freeze_sha256'],data_archive_sha256=data['archive_sha256'],retained_actual_pinned_comparator_replay_sha256=sha(d/'original_pinned_frozen_comparator_replay.json'),retained_pinned_compute_replay_verified=True,actual_frozen_comparator_executed_by_this_script=False,original_runtime_or_controller_guards_executed_by_this_script=False,private_guard_replacement_scope='Explicit public compressed source/data hash closure plus retained actual pinned-compute receipts. Only unchanged audit_math/audit_responses are called; main/closure/replay are never called.',resource_assessment=assessment,public_replay_script_sha256=sha(__file__),original_mathematical_source_changed=False,producer_eigensolver_sampler_or_new_fixture_executed=False,authorship=orig['authorship'])
  write(out/'public_replay_receipt.json',result);print(json.dumps(dict(status='PASS',role=role,public_checks=len(mathaudit.checks),closure_checks=len(checks),scientific_completions=result['scientific_completions'],receipt_sha256=sha(out/'public_replay_receipt.json')),indent=2))
if __name__=='__main__':main()
