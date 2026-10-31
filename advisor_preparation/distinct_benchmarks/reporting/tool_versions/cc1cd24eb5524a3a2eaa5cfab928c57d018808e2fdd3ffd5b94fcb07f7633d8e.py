#!/usr/bin/env python3
"""Audit immutable receipts and merge exactly 100 catalogue ideas; never run physics."""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math,re,shutil,tarfile,tomllib
from datetime import datetime,timezone
from pathlib import Path

STATUSES=('PASS','FAIL','UNRESOLVED','NOT_RUN')
ANALYTIC_FIXTURES=dict(zip(('A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 A11 A12 A13 A14 A22 A23').split(),[18,8,8,8,8,18,8,9,4,8,24,18,8,12,6,3]))
BASE=Path(__file__).resolve().parents[1]
REVIEWED_ED_PATCH_SHA256='da3fe6cde8c390f7dfb0dede91f8a04ea8e2afc04ae5c9448afe10d2b9c8c34f'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def finite(x):return isinstance(x,(float,int)) and not isinstance(x,bool) and math.isfinite(x)
def clean(x):
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,list):return [clean(v) for v in x]
 if isinstance(x,float) and not math.isfinite(x):return None
 return x
def dump(path,x):path.write_text(json.dumps(clean(x),indent=2,ensure_ascii=False)+'\n')
def metric_entry(name,error,tolerance,components=1):
 return {'metric':name,'maximum_error':error,'tolerance':tolerance,'components':components,'within_tolerance':finite(error) and finite(tolerance) and 0<=error<=tolerance}
def preserve_status(reported,failed=False,incomplete=False):
 if reported not in STATUSES:return 'UNRESOLVED'
 if reported!='PASS':return reported
 return 'FAIL' if failed else 'UNRESOLVED' if incomplete else 'PASS'
def source_scope(family,cid):
 if family=='journal':
  functions=['src/noninteracting.jl: bloch_hamiltonian']
  if cid=='J20':functions+=['sample_band_energies','zero_temperature_filling','src/realspace_kinetic.jl: realspace_kinetic_matrix','fourier_block']
  return {'method':'production free-band/real-space Hamiltonian builders','existing_user_functions':functions,'new_adapter':'published-equation evaluation, topology/contour or orbital-projector analysis','production_band_construction_exercised':True,'production_DQMC_sampling_exercised':False,'production_DMRG_exercised':False}
 if family=='packages':return {'method':'user ED state/Hamiltonian plus new observable adapter versus independent QuSpin','existing_user_functions':['analysis/compute_s2_ed.jl: build_space','apply_hamiltonian!','ground_state_lanczos']+(['renyi2_entropy_sites'] if cid=='P16' else []),'new_adapter':'explicit observable/partial-trace/spin-multiplet construction; see frozen package protocol','production_band_construction_exercised':False,'production_DQMC_sampling_exercised':False,'production_DMRG_exercised':False}
 ed=['build_space','apply_hamiltonian!'];thermal=['full_fock_matrices','exact_thermal_density']
 if cid in ('A03','A04','A05','A14'):method='ED Hamiltonian/Lanczos plus analytical adapter';functions=ed+['ground_state_lanczos']
 elif cid=='A22':method='ED Hamiltonian plus explicit excited-state/pair adapter';functions=ed
 elif cid=='A12':method='ED sector Hamiltonian and thermal ED density/RDM';functions=ed+['exact_thermal_density','reduced_operator']
 elif cid=='A23':method='thermal-ED Hamiltonian builder plus open-chain/full-spin adapter';functions=['full_fock_matrices']
 else:method='thermal ED plus analytical observable adapter';functions=thermal.copy()
 if cid in ('A05','A11'):functions+=['reduced_operator']
 if cid=='A05':functions+=['renyi2_entropy']
 if cid=='A13':functions+=['renyi2']
 return {'method':method,'existing_user_functions':functions,'existing_source_location':'analysis/compute_s2_ed.jl and/or analysis/audit_u4_l4_thermal_ed.jl','new_adapter':'sector, Lehmann, response, spectrum or theorem adapter declared in analytical protocol','production_band_construction_exercised':False,'production_DQMC_sampling_exercised':False,'production_DMRG_exercised':False}

def apply_single_file_patch(original,patch):
 """Apply checked unified hunks in memory; no source file is changed."""
 source=original.splitlines(keepends=True);lines=patch.splitlines(keepends=True);out=[];position=0;i=0;hunks=0
 while i<len(lines):
  line=lines[i]
  if not line.startswith('@@ '):i+=1;continue
  match=re.match(r'@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@',line)
  if not match:raise ValueError('Malformed unified patch hunk')
  start=int(match[1])-1;old_count=int(match[2] or 1);new_count=int(match[4] or 1)
  if start<position:raise ValueError('Overlapping/out-of-order patch hunks')
  out+=source[position:start];position=start;old_seen=new_seen=0;i+=1;hunks+=1
  while i<len(lines) and not lines[i].startswith('@@ '):
   line=lines[i]
   if line.startswith('\\ No newline'):i+=1;continue
   if not line or line[0] not in ' +-':raise ValueError('Unexpected unified patch line')
   if line[0] in ' -':
    if position>=len(source) or source[position]!=line[1:]:raise ValueError('Patch context differs from original source')
    position+=1;old_seen+=1
   if line[0] in ' +':out.append(line[1:]);new_seen+=1
   i+=1
  if (old_seen,new_seen)!=(old_count,new_count):raise ValueError('Unified patch hunk length differs')
 if not hunks:raise ValueError('No patch hunks')
 return ''.join(out+source[position:])

class Audit:
 def __init__(self,out):self.out=out;self.inputs={};self.warnings=[];self.parse_errors={}
 def input(self,role,path,required=False):
  path=Path(path).resolve()
  if not path.is_file():
   if required:raise ValueError(f'Required immutable input missing: {path}')
   self.warnings.append(f'Missing {role}: {path}');return None
  dest=self.out/'inputs'/f'{role}_{path.name}';dest.parent.mkdir(exist_ok=True)
  shutil.copyfile(path,dest)
  try:logical_path=path.relative_to(BASE.parent.resolve()).as_posix()
  except ValueError:logical_path=f'{role}/{path.name}'
  entry={'original_path':logical_path,'evidence_copy':str(dest.relative_to(self.out)),'bytes':path.stat().st_size,'sha256':sha(path)}
  if sha(dest)!=entry['sha256']:raise ValueError(f'Input changed while copying: {path}')
  self.inputs[role]=entry
  try:return json.loads(dest.read_text()) if path.suffix=='.json' else tomllib.loads(dest.read_text()) if path.suffix=='.toml' else dest
  except (ValueError,UnicodeError) as ex:
   self.parse_errors[role]=f'{type(ex).__name__}: {ex}';self.warnings.append(f'Malformed {role}; exact bytes preserved: {ex}')
   if required:raise ValueError(f'Required input cannot be parsed: {path}: {ex}') from ex
   return None

def audit_analytic(row):
 issues=[];metrics={};bad=False;count=0;fixtures=row.get('fixtures',[])
 if row.get('status')=='PASS' and len(fixtures)!=ANALYTIC_FIXTURES.get(row['id']):issues.append('Declared analytical fixture grid is incomplete or unexpected.')
 for fi,f in enumerate(fixtures,1):
  if not f.get('checks'):issues.append(f'Fixture {fi} has no comparisons.')
  for c in f.get('checks',[]):
   count+=1;m=c.get('metric');a=c.get('observed');b=c.get('analytical_reference');tol=c.get('tolerance')
   valid=m in ('maximum_absolute','L1') and isinstance(a,list) and isinstance(b,list) and bool(a) and len(a)==len(b) and all(finite(x) for x in a+b) and finite(tol) and tol>=0
   if not valid:issues.append(f'Invalid metric/value schema in fixture {fi}, check {count}.');continue
   errors=[abs(x-y) for x,y in zip(a,b)];error=sum(errors) if m=='L1' else max(errors)
   bad|=error>tol
   if c.get('status')!=('PASS' if error<=tol else 'FAIL'):issues.append(f'Reported check status disagrees with raw values at fixture {fi}, check {count}.')
   if not finite(c.get('error')) or abs(c['error']-error)>1e-12*max(1.,error):issues.append(f'Reported error disagrees with recomputation at fixture {fi}, check {count}.')
   current=metrics.setdefault(m,{'metric':m,'maximum_error':0.,'comparison_count':0,'scalar_components':0,'maximum_tolerance_ratio':0.,'all_within_tolerance':True})
   current['maximum_error']=max(current['maximum_error'],error);current['comparison_count']+=1;current['scalar_components']+=len(a)
   current['maximum_tolerance_ratio']=max(current['maximum_tolerance_ratio'],error/tol if tol else (0. if error==0 else math.inf));current['all_within_tolerance']&=error<=tol
 if row.get('status')=='PASS' and count==0:issues.append('PASS has no recorded comparisons.')
 return preserve_status(row.get('status'),bad,bool(issues)),list(metrics.values()),issues,len(fixtures),count

def audit_packages(row):
 issues=[];a=row.get('target_values');b=row.get('oracle_values');t=row.get('absolute_tolerances');bad=False;metrics=[]
 if row.get('status')=='PASS':
  valid=isinstance(a,list) and isinstance(b,list) and isinstance(t,list) and bool(a) and len(a)==len(b)==len(t) and all(finite(x) for x in a+b+t) and all(x>=0 for x in t)
  if not valid:issues.append('PASS lacks finite matched target/oracle values and component tolerances.')
  else:
   e=[abs(x-y) for x,y in zip(a,b)];bad=any(x>tol for x,tol in zip(e,t))
   metrics=[{'metric':'maximum_absolute','maximum_error':max(e),'scalar_components':len(e),'all_within_tolerance':not bad,'maximum_tolerance_ratio':max(x/tol if tol else (0. if x==0 else math.inf) for x,tol in zip(e,t))}]
   if not finite(row.get('maximum_absolute_error')) or abs(row['maximum_absolute_error']-max(e))>1e-12:issues.append('Reported package maximum disagrees with recomputation.')
 elif isinstance(a,list) and isinstance(b,list) and len(a)==len(b) and all(finite(x) for x in a+b):
  metrics=[{'metric':'maximum_absolute','maximum_error':max((abs(x-y) for x,y in zip(a,b)),default=0.),'scalar_components':len(a),'all_within_tolerance':False}]
 return preserve_status(row.get('status'),bad,bool(issues)),metrics,issues,1 if isinstance(a,list) else 0,len(a) if isinstance(a,list) else 0

def audit_journal(row):
 issues=[];metrics=[];bad=False;cid=row['id']
 def add(name,e,tol,n=1):
  nonlocal bad
  v=metric_entry(name,e,tol,n);metrics.append(v)
  if not finite(e):issues.append('Missing/nonfinite journal metric: '+name)
  else:bad|=not v['within_tolerance']
 if row.get('status') not in ('PASS','FAIL'):return row.get('status','UNRESOLVED'),metrics,issues,0,0
 if cid=='J20':
  meshes=row.get('meshes',[]);rays=row.get('fermi_rays',[])
  if [r.get('Nk') for r in meshes]!=[128,256,512,1024] or len(rays)!=64:issues.append('J20 registered mesh/ray grid incomplete.')
  for r in meshes:
   if r.get('published_topology_matches') is not True:bad=True
   add('band_energy_maximum_absolute',r.get('maximum_band_error'),1e-10)
   add('filling_maximum_absolute',r.get('filling_error'),1e-12)
  add('contour_radius_maximum_absolute',row.get('maximum_contour_radius_error'),1e-9,len(rays))
  if rays:
   actual=max(abs(x['target_radius']-x['published_equation_radius']) for x in rays)
   if abs(actual-row.get('maximum_contour_radius_error',math.inf))>1e-12:issues.append('J20 contour maximum disagrees with raw radii.')
  add('realspace_Fourier_block_maximum_absolute',row.get('realspace_Fourier_block_error'),1e-10)
  add('last_mesh_filling_change',row.get('last_mesh_filling_change'),1e-3)
  fixtures=len(meshes)+len(rays)
 elif cid=='J21':
  values=row.get('rows',[])
  if len(values)!=516:issues.append('J21 registered 2 paths x129 points x2 bands grid incomplete.')
  errors=[];sums=[]
  for r in values:
   a=r.get('target_orbital_weights',[]);b=r.get('published_projector_weights',[])
   if len(a)!=2 or len(b)!=2 or not all(finite(x) for x in a+b):issues.append('Invalid J21 orbital weights.');continue
   errors.append(max(abs(x-y) for x,y in zip(a,b)));sums.append(abs(sum(a)-1.))
  add('orbital_weight_maximum_absolute',max(errors,default=math.inf),1e-8,len(values)*2)
  add('weight_normalization_maximum_absolute',max(sums,default=math.inf),1e-12,len(values))
  if row.get('all_registered_orbital_transfers') is not True:bad=True
  if not finite(row.get('minimum_band_splitting')) or row['minimum_band_splitting']<=2e-8:issues.append('J21 nondegenerate-band guard is not established.')
  fixtures=len(values)
 else:issues.append('No frozen journal audit schema for this ID.');fixtures=0
 return preserve_status(row.get('status'),bad,bool(issues)),metrics,issues,fixtures,len(metrics)

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--result-dir',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True)
 ap.add_argument('--catalogue',type=Path,default=BASE.parent/'external_benchmark_roadmap'/'benchmark_ideas_100.json')
 ap.add_argument('--freeze',type=Path,default=BASE/'source_releases'/'tranche01-freeze.json')
 ap.add_argument('--source-archive',type=Path,default=BASE/'source_releases'/'tranche01-source.tar.gz')
 ap.add_argument('--journal',type=Path,help='Explicit selected later journal receipt; the original is preserved as history.')
 ap.add_argument('--journal-driver',type=Path,help='Exact selected journal driver; a v2 serialization-only amendment is verified against frozen v1.')
 ap.add_argument('--journal-predecessor',type=Path,action='append',default=[],help='Additional failed/earlier journal receipt; retained without contributing pass counts.')
 ap.add_argument('--packages',type=Path,help='Explicit selected later complete package-comparison receipt.')
 ap.add_argument('--package-ed-source',type=Path,help='Exact minimally patched ED source used by the selected complete seven-ID package suite.')
 ap.add_argument('--package-ed-patch',type=Path,help='Reviewed sparse P16 unified patch; applied only in memory for provenance verification.')
 ap.add_argument('--package-oracle-dir',type=Path,help='Immutable independent oracle directory; defaults to the original attempt oracle.')
 ap.add_argument('--package-predecessor',type=Path,action='append',default=[],help='Additional earlier package receipt; retained as history, never selected automatically.')
 ap.add_argument('--renyi-source-root',type=Path,default=BASE.parent/'validated_package'/'HubbardRenyiDQMC');a=ap.parse_args()
 out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False);audit=Audit(out)
 catalogue=audit.input('catalogue',a.catalogue,True);freeze=audit.input('frozen_release',a.freeze,True);archive=audit.input('frozen_archive',a.source_archive,True)
 ids=[c['id'] for c in catalogue]
 if len(ids)!=100 or len(set(ids))!=100:raise ValueError('Catalogue must contain exactly 100 unique ideas.')
 release_issues=[]
 if sha(archive)!=freeze.get('archive_sha256'):release_issues.append('Source archive SHA256 differs from frozen release.')
 archived_journal=None;analytic_grids={};package_manifest={}
 with tarfile.open(archive) as tar:
  for name,digest in freeze['files'].items():
   try:blob=tar.extractfile('distinct_benchmarks/'+name).read()
   except (KeyError,AttributeError):release_issues.append('Frozen source archive lacks '+name);continue
   if hashlib.sha256(blob).hexdigest()!=digest:release_issues.append('Archived source hash mismatch: '+name)
   if name=='journal/run_raghu.jl':archived_journal=blob
   if name=='analytic/manifest.toml':analytic_grids={x['id']:x['grid'] for x in tomllib.loads(blob.decode())['benchmarks']}
   if name=='packages/manifest.toml':package_manifest=tomllib.loads(blob.decode())
 expected=set(freeze['distinct_ids'])
 if not expected.issubset(ids) or len(expected)!=len(freeze['distinct_ids']):release_issues.append('Frozen distinct-ID list is invalid.')
 expected_renyi={}
 for role,name in [('user_ed_source','analysis/compute_s2_ed.jl'),('user_thermal_source','analysis/audit_u4_l4_thermal_ed.jl')]:
  source=audit.input(role,a.renyi_source_root/name,True);expected_renyi[name]=sha(source)
  if role=='user_ed_source':original_ed_copy=source
 package_ed_sha=expected_renyi['analysis/compute_s2_ed.jl'];package_amendment=None
 if a.package_ed_source:
  selected_source=audit.input('selected_package_ED_source',a.package_ed_source,True);package_ed_sha=sha(selected_source)
  if not a.package_ed_patch:raise ValueError('A selected patched package source requires --package-ed-patch')
  patch=audit.input('selected_package_ED_patch',a.package_ed_patch,True)
  transformed=apply_single_file_patch(original_ed_copy.read_text(),patch.read_text())
  package_amendment={'original_ED_sha256':expected_renyi['analysis/compute_s2_ed.jl'],'selected_ED_sha256':package_ed_sha,'patch_sha256':sha(patch),'reviewed_sparse_patch_verified':sha(patch)==REVIEWED_ED_PATCH_SHA256,'exact_source_transformation_verified':transformed==selected_source.read_text(),'scope':'Sparse fermionic subsystem regrouping repair; requires the complete unchanged package target/oracle/comparator suite, not an isolated candidate check.'}
 journal_driver_sha=freeze['files'].get('journal/run_raghu.jl');journal_amendment=None
 if a.journal_driver:
  driver=audit.input('selected_journal_driver',a.journal_driver,True);journal_driver_sha=sha(driver)
  if journal_driver_sha!=freeze['files'].get('journal/run_raghu.jl'):
   original=archived_journal.decode() if archived_journal else ''
   expected_text=original.replace('"center"=>center,"band"=>band','"center"=>Float64.(center),"band"=>band')
   journal_amendment={'original_driver_sha256':freeze['files'].get('journal/run_raghu.jl'),'selected_driver_sha256':journal_driver_sha,'change':'Recorded center arrays converted to Float64 for TOML serialization; no physics/gate change.','serialization_only_transformation_verified':bool(original) and expected_text!=original and driver.read_text()==expected_text}
 predecessors=[('journal',p) for p in a.journal_predecessor]+[('packages',p) for p in a.package_predecessor]
 if a.journal and a.journal.resolve()!=(a.result_dir/'journal.toml').resolve() and (a.result_dir/'journal.toml').exists():predecessors.insert(0,('journal',a.result_dir/'journal.toml'))
 if a.packages and a.packages.resolve()!=(a.result_dir/'packages.json').resolve() and (a.result_dir/'packages.json').exists():predecessors.append(('packages',a.result_dir/'packages.json'))
 history=[];history_paths=set()
 for family,path in predecessors:
  if path.resolve() in history_paths:continue
  history_paths.add(path.resolve());role=f'{family}_predecessor_{len(history)+1:02d}';previous=audit.input(role,path)
  previous_rows=previous.get('cases',previous.get('benchmarks',[])) if isinstance(previous,dict) else []
  history.append({'input_role':role,'family':family,'status':'MALFORMED' if role in audit.parse_errors else 'MISSING' if previous is None else 'VALID_NOT_SELECTED','reported_statuses':{r['id']:r.get('status') for r in previous_rows if isinstance(r,dict) and 'id' in r},'parse_error':audit.parse_errors.get(role),'counted_as_current_result':False})
 ledger={c['id']:{'id':c['id'],'category':c['category'],'title':c['title'],'scientific_question':c['scientific_question'],'catalogue_model':c['model_geometry'],'catalogue_ensemble':c['ensemble'],'source_urls':c['source_urls'],'status':'NOT_RUN','reported_status':None,'reason':'No execution receipt for this distinct idea. Earlier parameter-comparison batches do not establish completion.','fixtures':0,'comparisons':0,'metrics':[],'scope':{},'input_role':None,'audit_issues':[]} for c in catalogue}
 seen=set();raw_metadata={};all_family_issues={}
 for family,filename,key,validator in [('analytic','analytic.toml','results',audit_analytic),('packages','packages.json','cases',audit_packages),('journal','journal.toml','benchmarks',audit_journal)]:
  path=a.journal if family=='journal' and a.journal else a.packages if family=='packages' and a.packages else a.result_dir/filename
  raw=audit.input(family,path);family_issues=list(release_issues)
  if raw is None:continue
  if not isinstance(raw,dict):audit.warnings.append(f'{family} receipt is not an object; no completion inferred.');continue
  raw_metadata[family]={k:v for k,v in raw.items() if k!=key}
  if family=='analytic':
   for name,digest in raw.get('tranche_sha256',{}).items():
    if digest!=freeze['files'].get('analytic/'+name):family_issues.append('Analytical tranche source mismatch: '+name)
   if set(raw.get('tranche_sha256',{}))!={'run_analytic_tranche.jl','manifest.toml','PROTOCOL.md'}:family_issues.append('Incomplete analytical source closure.')
   if raw.get('source_sha256')!=expected_renyi:family_issues.append('Analytical imported user source differs from the audited source copy.')
  elif family=='packages':
   meta=raw.get('target_metadata',{});oracle=raw.get('oracle_receipt',{})
   for value,name in [(raw.get('manifest_sha256'),'packages/manifest.toml'),(raw.get('comparison_source_sha256'),'packages/compare.py'),(meta.get('adapter_sha256'),'packages/user_ed_distinct.jl'),(oracle.get('source_sha256'),'packages/quspin_distinct.py')]:
    if value!=freeze['files'].get(name):family_issues.append('Package provenance mismatch: '+name)
   if meta.get('user_solver_source_sha256')!=package_ed_sha:family_issues.append('Package target ED source mismatch.')
   if package_amendment and not (package_amendment['reviewed_sparse_patch_verified'] and package_amendment['exact_source_transformation_verified']):family_issues.append('Selected package source is not the independently reviewed exact sparse patch transformation.')
   for receipt_manifest in (meta.get('manifest_sha256'),oracle.get('manifest_sha256')):
    if receipt_manifest!=freeze['files'].get('packages/manifest.toml'):family_issues.append('Package target/oracle manifest differs from frozen protocol.')
   package_dir=Path(path).parent;oracle_dir=a.package_oracle_dir or a.result_dir/'oracle'
   for role,input_path,field in [('package_target',package_dir/'target.toml','target_file_sha256'),('package_oracle',oracle_dir/'oracle.toml','oracle_file_sha256')]:
    value=audit.input(role,input_path)
    if value is None or audit.inputs[role]['sha256']!=raw.get(field):family_issues.append('Package raw immutable input missing/mismatched: '+role)
   if not finite(meta.get('ground_energy_uncentered')) or not finite(oracle.get('ground_energy_uncentered')) or abs(meta['ground_energy_uncentered']-oracle['ground_energy_uncentered'])>1e-8:family_issues.append('Independent package ground-energy certification failed.')
  else:
   if raw.get('driver_sha256')!=journal_driver_sha:family_issues.append('Journal driver source mismatch.')
   if journal_amendment and not journal_amendment['serialization_only_transformation_verified']:family_issues.append('Selected journal source amendment is not the reviewed serialization-only transformation.')
   for name,digest in raw.get('source_hashes',{}).items():
    if digest!=freeze['files'].get('frozen_alm/'+name):family_issues.append('Production band source mismatch: '+name)
   if set(raw.get('source_hashes',{}))!={'src/noninteracting.jl','src/realspace_kinetic.jl'}:family_issues.append('Incomplete journal production-source closure.')
  all_family_issues[family]=family_issues
  for row in raw.get(key,[]):
   if not isinstance(row,dict):audit.warnings.append(f'Non-object {family} result excluded; no completion inferred.');continue
   cid=row.get('id')
   if cid not in ledger:audit.warnings.append(f'Unknown receipt ID {cid!r}; excluded from idea counts.');continue
   if cid in seen:
    ledger[cid]['status']='FAIL' if ledger[cid]['status']=='FAIL' or row.get('status')=='FAIL' else 'UNRESOLVED';ledger[cid]['audit_issues'].append('Duplicate receipt result; no favorable selection among duplicates is allowed.');continue
   seen.add(cid)
   try:status,metrics,issues,fixtures,comparisons=validator(row)
   except (ValueError,TypeError,KeyError,AttributeError) as ex:status='FAIL' if row.get('status')=='FAIL' else 'UNRESOLVED';metrics=[];issues=[f'Malformed result schema: {type(ex).__name__}: {ex}'];fixtures=comparisons=0
   issues+=family_issues
   if row.get('status')=='PASS' and cid not in expected:issues.append('PASS ID was not part of the frozen release.');status='UNRESOLVED'
   if issues and status=='PASS':status='UNRESOLVED'
   ledger[cid].update(status=status,reported_status=row.get('status'),reason=row.get('reason',''),fixtures=fixtures,comparisons=comparisons,metrics=metrics,scope=source_scope(family,cid),input_role=family,audit_issues=issues)
   if family=='analytic':ledger[cid]['executed_fixture_scope']=analytic_grids.get(cid)
   elif family=='packages':ledger[cid]['executed_fixture_scope']={k:package_manifest.get(k) for k in ('L','t','U','mu_uncentered','boundary','Nup','Ndown')}|{'observable':row.get('scientific_question'),'source_amendment_selected':bool(package_amendment)}
   else:ledger[cid]['executed_fixture_scope']={'scope':row.get('scope'),'parameters':raw.get('parameters'),'reference_version':raw.get('reference_version')}
 for cid in expected-seen:
  family='analytic' if cid.startswith('A') else 'journal' if cid.startswith('J') else 'packages'
  detail=audit.parse_errors.get(family,'Frozen attempted idea has no result in the provided receipts.')
  ledger[cid].update(status='UNRESOLVED',reason=detail+' Completion is not established.')
 process=audit.input('process_status',a.result_dir/'process_status.txt')
 if process is not None:raw_metadata['process_status_text']=process.read_text()
 rows=[ledger[cid] for cid in ids];counts={s:sum(r['status']==s for r in rows) for s in STATUSES}
 categories={cat:{s:sum(r['category']==cat and r['status']==s for r in rows) for s in STATUSES} for cat in sorted({r['category'] for r in rows})}
 summary={'generated_utc':datetime.now(timezone.utc).isoformat(),'tool_sha256':sha(__file__),'counts':counts,'category_counts':categories,'counting_rule':'Exactly 100 catalogue IDs; one status per distinct idea. Parameters, matrix elements, fixture/check/component counts are not additional ideas.','scope_limit':'ED/thermal-ED analysis and production free-band construction only for this tranche. No production DQMC sampling or DMRG pass is inferred.','release':freeze['release'],'release_archive_sha256':freeze['archive_sha256'],'release_audit_issues':release_issues,'journal_source_amendment':journal_amendment,'package_source_amendment':package_amendment,'predecessor_history':history,'parse_errors':audit.parse_errors,'input_manifest':audit.inputs,'raw_metadata':raw_metadata,'warnings':audit.warnings,'rows':rows}
 dump(out/'ledger_100.json',summary);dump(out/'input_manifest.json',audit.inputs)
 columns=['id','category','title','status','reported_status','input_role','fixtures','comparisons','metric_maxima','method_scope','reason','audit_issues']
 with (out/'ledger_100.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=columns);w.writeheader()
  for r in rows:w.writerow({k:r.get(k) for k in columns if k not in ('metric_maxima','method_scope','audit_issues')}|{'metric_maxima':json.dumps(r['metrics'],ensure_ascii=False),'method_scope':r['scope'].get('method',''),'audit_issues':' '.join(r['audit_issues'])})
 lines=['# Distinct scientific benchmark ledger','',f"**{counts['PASS']}/100 PASS; {counts['FAIL']} FAIL; {counts['UNRESOLVED']} UNRESOLVED; {counts['NOT_RUN']} NOT_RUN.**",'',summary['counting_rule'],'',summary['scope_limit'],'','Every raw input is copied byte-for-byte under `inputs/` and hashed in `input_manifest.json`. Inputs and the catalogue are never modified. PASS requires the frozen-source provenance and a completed scientific comparison. FAIL and incomplete/unresolved outcomes remain visible. A11 is summarized by its actual L1 metric, never the analytical driver’s empty generic maximum-absolute field.','','| Category | PASS | FAIL | UNRESOLVED | NOT_RUN |','|---|---:|---:|---:|---:|']
 for cat,n in categories.items():lines.append(f"| {cat} | {n['PASS']} | {n['FAIL']} | {n['UNRESOLVED']} | {n['NOT_RUN']} |")
 lines+=['','## Executed and unresolved ideas','','| ID | Status | Target scope | Metric-specific maxima |','|---|---|---|---|']
 for r in rows:
  if r['status']=='NOT_RUN':continue
  maxima='; '.join(f"{m['metric']}: {m['maximum_error']:.6g}" if finite(m.get('maximum_error')) else f"{m['metric']}: nonfinite/unavailable" for m in r['metrics']) or 'No completed numeric comparison'
  lines.append(f"| {r['id']} | {r['status']} | {r['scope'].get('method','Completion not established')} | {maxima} |")
 lines+=['','The complete 100-row ledger and all explanations are in [ledger_100.csv](ledger_100.csv) and [ledger_100.json](ledger_100.json). Fixture and scalar-component totals are descriptive evidence, not pass-count multipliers. Existing user analysis functions and new observable adapters are explicitly identified in the JSON; only J20/J21 exercise the original production band/real-space builders.','','## Provenance','',f"- Frozen release: `{freeze['release']}`",f"- Frozen archive SHA256: `{freeze['archive_sha256']}`",f"- Catalogue SHA256: `{audit.inputs['catalogue']['sha256']}`",f"- Reporting tool SHA256: `{summary['tool_sha256']}`"]
 issues=release_issues+audit.warnings+[f"{family}: {x}" for family,v in all_family_issues.items() for x in v]
 if issues:lines+=['','Audit notes:','']+['- '+x for x in issues]
 (out/'README.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'counts':counts,'output_dir':str(out),'warnings':len(issues)},indent=2))
if __name__=='__main__':main()
