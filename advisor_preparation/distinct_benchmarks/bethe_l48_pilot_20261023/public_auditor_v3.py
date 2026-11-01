#!/usr/bin/env python3
"""Public selection audit v3, derived from exact preserved v2; same quality gates."""
import argparse,hashlib,json,math,re,shutil,subprocess,sys,tarfile,tempfile,tomllib,traceback
from pathlib import Path
ARCHIVE='b1ebb51e6d6daa140e711fa08307d2384cceccdc78fcfa85db72d5ab778444ff'
FREEZE='4edfc75eff3e6b7b298e447999d72fbc6bddf0e70f2964242cf8e2129affe185'
PREDECESSOR='8c93bb6fe7f1128cc168f79b56d8b4366a7773058605aa0c3f2bba4b23914291'
JOB='6466266'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return tomllib.loads(p.read_text())if p.suffix=='.toml'else json.loads(p.read_text())
def finite(v):
 if isinstance(v,dict):return all(finite(x)for x in v.values())
 if isinstance(v,list):return all(finite(x)for x in v)
 return not isinstance(v,float)or math.isfinite(v)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--public-inputs',action='store_true');ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();raw=a.raw;out=a.output;assert not out.exists();checks=[];report={};rb=None;quality_failures=[]
 def require(b,label):checks.append({'name':label,'passed':bool(b),'kind':'integrity'});assert b,label
 def quality(b,label):
  checks.append({'name':label,'passed':bool(b),'kind':'quality'})
  if not b:quality_failures.append(label)
 try:
  source=raw/'source';cfg=load(source/'pilot.toml');fr=load(raw/'release_freeze.json');archive=raw/'source_archive.tar.gz';target=load(raw/'results/target.toml');summary=load(raw/'results/summary.json');oracle=load(raw/'results/oracle.json')
  require(sha(raw/'release_freeze.json')==FREEZE and sha(archive)==fr['archive_sha256']==ARCHIVE,'Exact original L48 source freeze/archive')
  with tarfile.open(archive)as tf:
   ms=tf.getmembers();require(len(ms)==len(fr['files'])==12,'Full12-member original source closure')
   for m in ms:
    require(m.isfile()and not Path(m.name).is_absolute()and'..'not in Path(m.name).parts,'Safe source archive member');n=str(Path(m.name).relative_to(fr['prefix'].rstrip('/')));b=tf.extractfile(m).read();require(hashlib.sha256(b).hexdigest()==fr['files'][n]and b==(source/n).read_bytes(),'Frozen archive/source identity '+n)
  require(cfg['L']==[48]and cfg['phase']=='PILOT_ONLY'and not cfg['scientific_completion_allowed'],'L48-only excluded pilot')
  require(sha(raw/'predecessor_auditor.py')==PREDECESSOR,'Preserved reviewed L32 auditor predecessor');priorcfg=load(raw/'predecessor_manifest.toml');require(sha(raw/'predecessor_manifest.toml')==cfg['predecessor_manifest_sha256'],'Exact prior L32 manifest');require(all(cfg[k]==priorcfg[k]for k in ['gates','arms','oracle','U','t','mu_physical','boundary','source_sha256','project_sha256','manifest_sha256','julia_version','itensors_version','itensormps_version','numpy_version','scipy_version']),'All unchanged prior model/precision/oracle/runtime gates')
  if a.public_inputs:
   inventory=load(raw/'public_input_inventory.json');require(inventory['schema']=='L48_PUBLIC_SCIENTIFIC_INPUT_SELECTION_V1' and inventory['job_id']==JOB,'Disclosed public scientific selection')
   require(inventory['original_audit_sha256']=='3a5a36ecc454e2a9a8f0118fb3394e061c20de1228b0c1ecec116021a46e03e0','Preserved complete private audit identity')
   for n,r in inventory['files'].items():
    p=raw/n;require(not Path(n).is_absolute() and '..' not in Path(n).parts and p.is_file() and not p.is_symlink() and p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],'Published scientific bytes/hash '+n)
   for role,r in inventory['postexecution_runtime_hashes'].items():require(len(r['sha256'])==64 and r['bytes']>0,'Disclosed postexecution runtime hash '+role)
   original_transfer_file_count=inventory['original_transfer_file_count']
  else:
   before=load(raw/'transfer_before.json');after=load(raw/'transfer_after.json');aft={r['local']:r for r in after['files']};require(before['job_id']==JOB and after['status']=='PASS','Original job and remote-after transfer')
   for r in before['files']:
    p=raw/r['local'];v=aft[r['local']];require(p.is_file()and not p.is_symlink()and p.stat().st_size==r['bytes']and sha(p)==r['sha256'],'Transferred bytes/hash '+r['local']);require(v['sha256']==r['sha256']and v['bytes']==r['bytes']and v['mtime_ns']==r['before']['mtime_ns']and v['inode']==r['before']['inode'],'Remote unchanged '+r['local'])
   for role,r in before['runtime_files'].items():require(len(r['sha256'])==64 and r['bytes']>0 and after['runtime_files'][role]['sha256']==r['sha256'],'Read-only post-run runtime hash inventory '+role)
   original_transfer_file_count=len(before['files'])
  require(finite(cfg)and finite(target)and finite(summary)and finite(oracle),'Finite complete numeric receipts')
  require(target['manifest_sha256']==oracle['manifest_sha256']==summary['manifest_sha256']==sha(source/'pilot.toml'),'Manifest bindings')
  require(target['adapter_sha256']==sha(source/'pilot_target.jl')and oracle['oracle_source_sha256']==sha(source/'oracle_lieb_wu.py')and summary['analysis_source_sha256']==sha(source/'summarize_pilot.py'),'Actual target/oracle/summarizer source hashes')
  for key,field,rel in [('source_sha256','source_sha256','analysis/compute_s2_itensor_dmrg.jl'),('project_sha256','project_sha256','Project.toml'),('manifest_sha256','environment_manifest_sha256','Manifest.toml')]:require(sha(raw/'native'/rel)==cfg[key]==target[field],'Actual user MPO/environment '+rel)
  require(target['julia_version']==cfg['julia_version']and target['ITensors_version']==cfg['itensors_version']and target['ITensorMPS_version']==cfg['itensormps_version'],'Native actual runtime versions')
  require(oracle['versions']=={k:cfg[k+'_version']for k in ['numpy','scipy']},'Independent numerical runtime versions')
  ep=cfg['engineering_prerequisite'];prior=source/ep['audit_relative'];anchor=source/ep['anchor_relative'];previous=load(prior);r=previous['resources'];ab=anchor.read_bytes();require(sha(prior)==ep['audit_sha256']and sha(anchor)==ep['preexecution_engineering_anchor_sha256'],'Original engineering prerequisite hashes')
  require(previous['status']=='PASS'and previous['checks']==177 and previous['scientific_completions']==0 and previous['auditor_sha256']==PREDECESSOR,'Preserved177-check finiteL32 audit')
  require(ep['preexecution_engineering_anchor_commit']=='66950c791f7f939ca176e5427be62cd8dfca8348'and hashlib.sha1(b'blob '+str(len(ab)).encode()+b'\0'+ab).hexdigest()=='6f09e3c24b16752df0bf6577ae9d9c631dc20af7','Preserved published engineering anchor identity')
  require('at most 30 minutes and peak RSS at most 2300 MiB'in ab.decode(),'Original published engineering criterion unchanged')
  require(r['job_elapsed_seconds']==ep['L32_job_elapsed_seconds']==1316<=ep['maximum_L32_wall_seconds']==1800,'Actual prior L32 elapsed/criterion')
  require(r['batch_peak_rss_kib']==ep['L32_batch_peak_rss_kib']==1735276 and r['target_process_peak_rss_kib']==ep['L32_target_peak_rss_kib']==1740640 and max(r['batch_peak_rss_kib'],r['target_process_peak_rss_kib'])<=ep['maximum_L32_peak_rss_MiB']*1024==2300*1024,'Actual prior L32 memory/criterion')
  require(target['status']=='PILOT_ONLY_COMPLETE'and not target['scientific_completion_allowed'],'Complete finiteL48 pilot, no thermodynamic claim')
  rows={(x['sector'],x['arm']):x for x in target['rows']};require(len(rows)==len(target['rows'])==6 and set(rows)=={(s,b)for s in cfg['sectors']for b in cfg['arms']},'All six sector/precision solves')
  g=cfg['gates']
  for(sector,arm),x in rows.items():
   require(load(raw/'results'/f'L48_{sector}_{arm}.toml')==x,'Individual/aggregate exact identity '+sector+'/'+arm);ac=cfg['arms'][arm];nu,nd=24+int(sector=='plus'),24-int(sector=='minus');require((x['L'],x['Nup_expected'],x['Ndown_expected'])==(48,nu,nd),'Expected fixed particle/spin sector '+sector+'/'+arm)
   require(x['seed']==ac['seed_base']+48000+10*(cfg['sectors'].index(sector)+1),'Original seed algorithm '+sector+'/'+arm);quality(abs(x['norm']-1)<=g['norm_abs'],'Norm '+sector+'/'+arm);quality(abs(x['Nup']-nu)<=g['sector_abs']and abs(x['Ndown']-nd)<=g['sector_abs'],'Actual particle/spin projection '+sector+'/'+arm)
   quality(abs(x['energy']-x['reported_energy'])<=g['reported_contracted_energy_abs'],'Reported/contracted physical energy '+sector+'/'+arm);quality(g['variance_min']<=x['variance']and x['variance']/48<=ac['variance_per_site_max'],'Original signed variance gates '+sector+'/'+arm);require(abs(x['variance_per_site']-x['variance']/48)<=1e-12 and abs(x['energy_per_site']-x['energy']/48)<=1e-12,'Per-site units '+sector+'/'+arm)
   require(len(x['sweep_truncerrors'])==ac['nsweeps']and x['final_truncerr']==x['sweep_truncerrors'][-1],'Complete unchanged sweep/truncation records '+sector+'/'+arm);quality(0<=x['final_truncerr']<=ac['final_truncerr_max'],'Original final truncation gate '+sector+'/'+arm);quality(0<x['max_bond']<=ac['chi_max'],'Declared bond cap '+sector+'/'+arm)
   if sector=='half':quality(g['double_occupancy_min']<=x['double_occupancy_per_site']<=g['double_occupancy_max'],'Half double-occupancy range '+arm)
   for name,total in [('site_Nup',x['Nup']),('site_Ndown',x['Ndown']),('site_double_occupancy',48*x['double_occupancy_per_site'])]:require(len(x[name])==48 and abs(sum(x[name])-total)<=1e-9,'All48 local values/normalization '+name+'/'+sector+'/'+arm)
  values={}
  for arm in cfg['arms']:
   h,m,q=(rows[s,arm]for s in ['half','minus','plus']);quality(abs(q['energy']-m['energy']-cfg['U'])<=g['particle_hole_energy_abs'],'Particle-hole Nplus/Nminus identity '+arm);values[arm]={'energy_per_site':h['energy']/48,'charge_gap':q['energy']+m['energy']-2*h['energy'],'double_occupancy':sum(h['site_double_occupancy'])/48}
  errors={k:abs(values['coarse'][k]-values['refined'][k])for k in values['coarse']}
  for key,v in errors.items():quality(v<=g['coarse_refined_'+key+'_abs'],'Unchanged coarse/refined stability '+key)
  require(oracle['status']=='COMPLETED'and oracle['oracle_gate']=='PASS','Original thermodynamic numerical reference eligible')
  for name,x in oracle['values'].items():
   require(x['estimated_absolute_error']<=cfg['oracle']['maximum_estimated_absolute_error'],'Numerical reference error budget '+name);require(all(abs(x['adaptive_integral']-z)<=cfg['oracle']['maximum_estimated_absolute_error']for z in x['fixed_gauss'].values()),'Adaptive/Gauss receipt agreement '+name)
  require(oracle['values']['energy']['value']==oracle['values']['energy']['adaptive_integral']and oracle['values']['double_occupancy']['value']==oracle['values']['double_occupancy']['adaptive_integral']and abs(oracle['values']['charge_gap']['value']-(cfg['U']-4*cfg['t']+oracle['values']['gap_integral']['adaptive_integral']))<=1e-14,'Physical reference conventions')
  require(abs(oracle['energy_derivative_five_point']-oracle['values']['double_occupancy']['value'])==oracle['energy_derivative_error']and oracle['energy_derivative_error']<=cfg['oracle']['energy_derivative_check_abs'],'Preserved derivative/D receipt identity')
  with tempfile.TemporaryDirectory(prefix='bethe_l48_postprocess_')as tmp:
   dest=Path(tmp)/'replay.json';subprocess.run([sys.executable,'-B',str(source/'summarize_pilot.py'),'--manifest',str(source/'pilot.toml'),'--oracle',str(raw/'results/oracle.json'),'--target',str(raw/'results/target.toml'),'--output',str(dest)],check=True,capture_output=True,text=True);rb=dest.read_bytes();require(rb==(raw/'results/summary.json').read_bytes(),'Byte-identical unchanged frozen summarizer replay')
  require(summary['pilot_quality']==('PILOT_QUALITY_UNRESOLVED' if quality_failures else 'PILOT_QUALITY_PASS')and summary['solves']==6 and summary['distinct_ideas_passed']==0 and {x['id']for x in summary['rows']}==set(cfg['catalogue_ids'])and all(x['status']=='UNRESOLVED'for x in summary['rows']),'Pilot quality and allthermodynamic IDs unexecuted')
  text=(raw/'results/target.resources.txt').read_text();rss=int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',text)[1]);require('Exit status: 0'in text and rss<=cfg['resources']['memory_mb']*1024 and target['elapsed_seconds']<cfg['resources']['target_command_cap_seconds'],'Bounded successful native process')
  accounting=(raw/'sacct_fresh.txt').read_text();jobs=[x.split('|')for x in accounting.splitlines()if x.startswith(JOB+'|')];require(len(jobs)==1,'One original scheduler job record');j=jobs[0];require(j[2:4]==['COMPLETED','0:0']and j[5:8]==['1','3500M','01:00:00'],'Original scheduler completion/cpu/memory/allocation');batch=[x.split('|')for x in accounting.splitlines()if x.startswith(JOB+'.batch|')];batchrss=batch[0][-1]if batch else None
  report.update(status='UNRESOLVED' if quality_failures else 'PASS',audit_integrity='PASS',pilot_quality=summary['pilot_quality'],quality_failures=quality_failures,original_quality_reasons=summary['quality_reasons'],refined=values['refined'],coarse_refined_errors=errors,resources={'job_elapsed_seconds':int(j[4]),'target_adapter_seconds':target['elapsed_seconds'],'target_process_peak_rss_KiB':rss,'batch_MaxRSS_original':batchrss,'frozen_source_wall_minutes':60,'actual_scheduler_wall_minutes':60,'target_cap_seconds':3300,'no_new_size_admission':True},source_archive_sha256=ARCHIVE,source_freeze_sha256=FREEZE,source_members=12,raw_text_files=original_transfer_file_count,native_source_environment_pins=3,runtime_hash_scope='Observed read-only postexecution runtime inventory plus source-checked native Project/Manifest and recorded actual versions; no retrospective preexecution binary-pin claim.')
 except Exception as ex:report.update(status='UNRESOLVED',audit_integrity='FAIL',quality_failures=quality_failures,reason=type(ex).__name__+': '+str(ex),traceback=traceback.format_exc())
 out.mkdir();shutil.copyfile(__file__,out/'auditor_source.py')
 if rb is not None:(out/'frozen_summary_replay.json').write_bytes(rb)
 inputs={str(f.relative_to(raw)):{'sha256':sha(f),'bytes':f.stat().st_size}for f in sorted(raw.rglob('*'))if f.is_file()and f.name!='remote_text_receipts.tar'}
 report.update(schema='BETHE_L48_PUBLIC_SELECTION_AUDIT_V3',public_input_selection=a.public_inputs,private_transport_checks_replayed=not a.public_inputs,checks=checks,check_count=len(checks),quality_failure_count=len(quality_failures),scientific_completions=0,catalogue_status={k:'NOT_RUN'for k in ['A18','A19','A20']},input_sha256={n:v['sha256']for n,v in inputs.items()},inputs=inputs,auditor_sha256=sha(__file__),adapted_from_preserved_L32_auditor_sha256=PREDECESSOR,physics_executed=False,oracle_or_DMRG_rerun=False,future_growth_permission=False,scope='Postprocessing only. Allquality/resource checks concern only this finiteL48 pilot; no thermodynamic equality/extrapolation/heldout validation and noL56/L64 admission.')
 (out/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');fs=[f for f in out.iterdir()if f.is_file()];(out/'SHA256SUMS').write_text(''.join(sha(f)+'  '+f.name+'\n'for f in sorted(fs)));print(json.dumps({k:report[k]for k in ['status','audit_integrity','check_count','quality_failure_count','scientific_completions']}))
if __name__=='__main__':main()
