#!/usr/bin/env python3
"""Saved-data J23 pilot audit. No HF solver, quadrature or scientific job is run."""
import argparse,hashlib,json,math,re,shutil,tarfile,tomllib,traceback
from pathlib import Path
ARCHIVE='5b78ce7ea96a31aedb5ee2c2a27e3f128667ae0121cf410219a13a4ae6ae6bbe'
FREEZE='7b763ad944dd745d604957bf25556fdad06b2ff733f09cb6c25230f7db7b4043'
JOB='6466377'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return tomllib.loads(p.read_text()) if p.suffix=='.toml' else json.loads(p.read_text())
def flatten(x):
 if isinstance(x,list):
  for v in x:yield from flatten(v)
 else:yield x

def independent_native(row,c):
 """Closed2x2 spectral projectors on the already recorded native momentum grid."""
 N=row['mesh'];m=row['m'];t4=row['t4'];U=2*m;align=row['alignment'];rho=[[0.,0.],[0.,0.]];kin=0.;bloch=[];eigs=[[],[]];occ=[[],[]];gap=math.inf
 V=[[-m,m],[-m,m]] if align=='parallel' else [[-m,m],[m,-m]]
 for i in range(N):
  kx=-math.pi+2*math.pi*(i+c['native_mesh_shift'])/N
  for j in range(N):
   ky=-math.pi+2*math.pi*(j+c['native_mesh_shift'])/N
   ex=-2*c['t1']*math.cos(kx)-2*c['t2']*math.cos(ky)-4*c['t3']*math.cos(kx)*math.cos(ky)
   ey=-2*c['t2']*math.cos(kx)-2*c['t1']*math.cos(ky)-4*c['t3']*math.cos(kx)*math.cos(ky)
   h=-4*t4*math.sin(kx)*math.sin(ky);bloch.append([ex,ey,h,0.])
   for spin in range(2):
    a=(ex+V[0][spin]+ey+V[1][spin])/2;z=(ex+V[0][spin]-ey-V[1][spin])/2;R=math.hypot(z,h)
    if not R>0:raise ValueError('Pilot unexpectedly has degenerate2x2 auxiliary block')
    energies=[a-R,a+R];f=[int(en<row['native_evaluate']['mu']) for en in energies];eigs[spin].append(energies);occ[spin].append(f)
    for sign,en,occupied in zip([-1,1],energies,f):
     gap=min(gap,abs(en-row['native_evaluate']['mu']))
     if occupied:
      p0=(1+sign*z/R)/2;p1=(1-sign*z/R)/2;p01=sign*h/(2*R)
      rho[0][spin]+=p0/(N*N);rho[1][spin]+=p1/(N*N);kin+=(ex*p0+ey*p1+2*h*p01)/(N*N)
 F=kin+U*sum((r[0]-.5)*(r[1]-.5) for r in rho)
 return dict(rho=rho,kinetic=kin,physicalF=F,V=V,bloch=bloch,eigenvalues=eigs,occupations=occ,gap=gap)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not a.output.exists();raw=a.raw;source=raw/'source';data=raw/'results';checks=[];metrics={};report={}
 def need(v,label):checks.append(dict(name=label,passed=bool(v)));assert v,label
 def close(x,y,lim,label):
  xs=list(flatten(x)) if isinstance(x,list) else [x];ys=list(flatten(y)) if isinstance(y,list) else [y]
  need(len(xs)==len(ys) and all(isinstance(v,(int,float)) and math.isfinite(v) for v in xs+ys),label+' finite/count')
  err=max((abs(x-y) for x,y in zip(xs,ys)),default=0.);need(err<=lim,label);metrics[label]=err
 try:
  fr=load(raw/'release_freeze.json');c=load(source/'manifest.toml');ex=load(data/'execution.json');comparison=load(data/'comparison.json');reference=load(data/'reference.json');fresh=load(a.replay/'comparison_replay.json');assessment=load(a.replay/'assessment.json')
  need(sha(raw/'release_freeze.json')==FREEZE and sha(raw/'source_archive.tar.gz')==fr['archive_sha256']==ARCHIVE,'Exact original frozen source/archive')
  with tarfile.open(raw/'source_archive.tar.gz') as tf:
   ms=tf.getmembers();need(len(ms)==len(fr['files'])==19,'All19 original source members')
   for m in ms:
    need(m.isfile() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts,'Safe archive member')
    n=str(Path(m.name).relative_to(fr['prefix'].rstrip('/')));b=tf.extractfile(m).read();need(hashlib.sha256(b).hexdigest()==fr['files'][n] and b==(source/n).read_bytes(),'Original frozen bytes '+n)
  need(sha(source/'source_freeze.json')==FREEZE,'Source-local original freeze')
  sourcehash={n:h for h,n in (line.split('  ',1) for line in (source/'SHA256SUMS').read_text().splitlines())}
  need(len(sourcehash)==18 and all(sha(source/n)==h for n,h in sourcehash.items()),'Exact18 substantive source hashes')
  before=load(raw/'transfer_before.json');after=load(raw/'transfer_after.json');aft={x['local']:x for x in after['files']}
  need(before['job_id']==JOB and after['status']=='PASS','Original pilot and unchanged remote-after receipts')
  for r in before['files']:
   p=raw/r['local'];q=aft[r['local']];need(p.is_file() and not p.is_symlink() and p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],'Transferred original bytes '+r['local']);need(q['bytes']==r['bytes'] and q['sha256']==r['sha256'] and q['mtime_ns']==r['before']['mtime_ns'] and q['inode']==r['before']['inode'],'Remote before/after unchanged '+r['local'])
  library=load(source/'reference_library_pins.json')
  for name,h in library['files'].items():
   role='reference/'+name;r=before['runtime_files'][role];q=after['runtime_files'][role];need(r['sha256']==q['sha256']==h and r['bytes']==q['bytes']==library['bytes'][name] and r['mtime_ns']==q['mtime_ns'] and r['inode']==q['inode'],'Actual39-file reference before/after pin '+name)
  jobrows=[line.split('|') for line in (raw/'sacct_fresh.txt').read_text().splitlines() if line.startswith(JOB+'|')];need(len(jobrows)==1 and jobrows[0][2:4]==['COMPLETED','0:0'] and jobrows[0][5:8]==['1','3500M','00:15:00'],'Original scheduler completion/CPU/memory/allocation')
  need(ex['status']=='PROCESSES_COMPLETE' and ex['mode']=='pilot' and ex['job_id']==JOB and ex['source_freeze_sha256']==FREEZE and ex['source_sha256']==sourcehash and ex['scientific_completions']==0 and not ex['SCF_executed'],'Original completed excluded pilot execution')
  names=['native16','native32','reference','comparator'];need([p['name'] for p in ex['processes']]==names and all(p['exit_code']==0 and not p['timed_out'] and math.isfinite(p['elapsed_seconds']) and p['elapsed_seconds']>0 for p in ex['processes']),'All four original processes completed')
  need(ex['cpus']==1 and ex['memory_MB']==3500 and ex['finished_unix']-ex['started_unix']<c['resources']['pilot_combined_timeout_seconds'],'Original bounded oneCPU pilot')
  for n,r in ex['outputs'].items():need((data/n).stat().st_size==r['bytes'] and sha(data/n)==r['sha256'],'Execution-bound scientific output '+n)
  need(reference['schema']=='J23_INDEPENDENT_ANALYTIC_WICK_RAW_V1' and reference['mode']=='pilot' and reference['job_id']==JOB and reference['source_sha256']==sourcehash and reference['source_freeze_sha256']==FREEZE and reference['reference_library_pins_sha256']==sha(source/'reference_library_pins.json') and reference['reference_library_sha256']==load(source/'reference_library_pins.json')['files'],'Original reference source/runtime pins')
  need(reference['python_version']==c['python_version'] and reference['versions']=={'numpy':c['numpy_version'],'scipy':c['scipy_version']} and not reference['SCF_executed'] and not reference['ground_state_claim'],'Pinned reference version and trial-only scope')
  allrows=[];lim=c['gates']
  for N in c['pilot_meshes']:
   doc=load(data/f'native{N}.toml');need(doc['mode']=='pilot' and doc['schema']=='J23_ACTUAL_NATIVE_HF_RAW_V1' and doc['job_id']==JOB and doc['source_sha256']==sourcehash and doc['source_freeze_sha256']==FREEZE and doc['native_source_sha256']==c['native_source_hashes'] and doc['julia_version']==c['julia_version'] and doc['trial_state_only'] and not doc['SCF_executed'],'Actual native original source/version '+str(N))
   need(doc['Ns']==1 and doc['d']==2 and len(doc['grids'])==1 and doc['grids'][0]['mesh']==N and doc['grids'][0]['t4']==c['pilot_t4'],'Actual physical-cell/grid normalization '+str(N))
   expected={(phase,align) for phase in c['pilot_phases'] for align in c['alignments']};need(len(doc['rows'])==4 and {(r['phase'],r['alignment']) for r in doc['rows']}==expected,'Complete four cold/warm rows '+str(N))
   for row in doc['rows']:
    tag=f"N{N}/{row['phase']}/{row['alignment']}";need(row['mesh']==N and row['m']==c['pilot_m'] and row['t4']==c['pilot_t4'] and row['U']==2*c['pilot_m'] and row['Ns']==1 and row['d']==2 and row['nk']==N*N,'Fixed row '+tag);own=independent_native(row,c);actual=row['native_evaluate'];obs=row['observer']
    close(obs['bloch_entries_by_k'],own['bloch'],1e-10,'Independent complete Bloch '+tag);close(obs['eigenvalues_by_spin_k_band'],own['eigenvalues'],1e-10,'Independent complete eigenvalues '+tag);close(obs['occupations_by_spin_k_band'],own['occupations'],0.,'Independent complete pure occupations '+tag);close(actual['V'],own['V'],1e-10,'Independent exchange spin signs '+tag)
    need(own['gap']>=lim['common_mu_gap_min'] and obs['pure_exact_binary_occupations'] and obs['two_occupied_per_k'],'Pure gapped trial qualification '+tag);close(obs['occupations_per_k'],[2.]*(N*N),0.,'All momentum counts '+tag);close(actual['n'],2.,lim['native_density_abs'],'Actual density '+tag);close(actual['entropy'],0.,lim['native_entropy_abs'],'Actual entropy '+tag)
    close(actual['rho'],own['rho'],lim['observer_density_vs_actual_native_rho_abs'],'Independent actual rho '+tag);close(obs['rho'],actual['rho'],lim['observer_density_vs_actual_native_rho_abs'],'Observer actual rho '+tag)
    for field in ['physicalF','kinetic']:
     close(actual[field],own[field],lim['same_mesh_physical_energy_per_cell_abs'],'Independent physical Wick '+field+'/'+tag);close(obs[field],actual[field],lim['same_mesh_physical_energy_per_cell_abs'],'Observer actual '+field+'/'+tag)
    need(math.isfinite(obs['maximum_eigenpair_residual']) and 0<=obs['maximum_eigenpair_residual']<=1e-10,'Eigenpair residual '+tag);allrows.append((row,own))
  refs={(r['grid_kind'],r['size'],r['alignment']):r for r in reference['rows']};expected={(kind,N,align) for kind,ns in [('native_mesh',c['pilot_meshes']),('gauss',c['oracle_gauss_orders'])] for N in ns for align in c['alignments']};need(len(reference['rows'])==len(refs)==10 and set(refs)==expected,'All four native-reference and six Gaussian rows')
  for k,r in refs.items():
   tag=str(k);need(r['m']==c['pilot_m'] and r['t4']==c['pilot_t4'] and r['U']==2*c['pilot_m'] and r['nk']==r['size']**2 and r['pure_occupations_qualified'] and r['common_mu_gap_min']>=1.,'Reference qualification '+tag);close(r['quadrature_weight_sum'],1.,1e-12,'Reference measure '+tag);close(r['density'],2.,1e-12,'Reference density '+tag);need(all(0<=v<=1 for v in flatten(r['rho'])),'Physical local reference density '+tag)
   wick=r['kinetic']+r['U']*sum((x[0]-.5)*(x[1]-.5) for x in r['rho']);close(r['physicalF'],wick,1e-12,'Reference physical Wick identity '+tag)
   if r['alignment']=='parallel':close(r['rho'],[[1.,0.],[1.,0.]],1e-12,'Parallel filled-spin projector '+tag);close(r['kinetic'],0.,1e-12,'Parallel full-band kinetic trace '+tag);close(r['physicalF'],-c['pilot_m'],1e-12,'Parallel analytic physical energy '+tag)
  for row,own in allrows:
   r=refs['native_mesh',row['mesh'],row['alignment']];tag=f"{row['mesh']}/{row['phase']}/{row['alignment']}"
   close(r['rho'],own['rho'],1e-10,'Independent native-reference rho '+tag);close(r['kinetic'],own['kinetic'],1e-10,'Independent native-reference kinetic '+tag);close(r['physicalF'],own['physicalF'],1e-10,'Independent native-reference physicalF '+tag)
  need(comparison['input_sha256']=={n:sha(data/n) for n in ['native16.toml','native32.toml','reference.json']},'Comparator original raw input hashes')
  need(comparison['status']=='PASS' and comparison['mode']=='pilot' and comparison['catalogue_status']=='NOT_RUN' and comparison['scientific_completions']==0 and comparison['pilot_excluded'] and comparison['source_freeze_sha256']==FREEZE and all(r['passed'] for r in comparison['checks']) and comparison['native_rows']==8 and comparison['reference_rows']==10,'Original full excluded-pilot comparator gates')
  # Fixed replay policy admits roundoff only in numerical gate values; decisions/bounds/status remain exact.
  def eq(x,y,path=()):
   if isinstance(x,dict):return isinstance(y,dict) and x.keys()==y.keys() and all(eq(v,y[k],path+(k,)) for k,v in x.items())
   if isinstance(x,list):return isinstance(y,list) and len(x)==len(y) and all(eq(v,w,path+(i,)) for i,(v,w) in enumerate(zip(x,y)))
   if isinstance(x,float) and path[-1:]==('value',) and 'checks' in path:return isinstance(y,float) and math.isfinite(x) and math.isfinite(y) and abs(x-y)<=max(1e-12,64*2.220446049250313e-16*max(1.,abs(x),abs(y)))
   return type(x)==type(y) and x==y
  need(eq(comparison,fresh),'Unchanged frozen comparator replay identities/decisions and declared arithmetic tolerance')
  need(assessment['schema']=='J23_EXCLUDED_PILOT_ENGINEERING_V1' and assessment['source_freeze_sha256']==FREEZE and assessment['execution_sha256']==sha(data/'execution.json') and assessment['comparison_sha256']==sha(data/'comparison.json') and assessment['pilot_job_id']==JOB and assessment['scientific_completions']==0 and not assessment['SCF_executed'],'Unchanged frozen assessor original inputs and zero credit')
  for name,h in assessment['input_sha256'].items():need((data/name).is_file() and sha(data/name)==h,'Frozen assessor original input '+name)
  times={p['name']:p['elapsed_seconds'] for p in ex['processes']};rss=[]
  for name in names:
   text=(data/(name+'.resources.txt')).read_text();hit=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',text);need(hit is not None and 'Exit status: 0' in text,'Original GNUtime success '+name);rss.append(int(hit[1]))
  n32=load(data/'native32.toml');warm=sum(r['elapsed_seconds'] for r in n32['rows'] if r['phase']=='warm')+n32['grids'][0]['grid_seconds'];pairs=len(c['m_values'])*len(c['t4_values']);sumN=sum(N*N for N in c['native_square_mesh']);pilotN=2*sum(N*N for N in c['pilot_meshes']);fullN=2*pairs*sumN;gaussN=sum(N*N for N in c['oracle_gauss_orders'])
  components=dict(native=2*pairs*sumN/(32*32)*warm+2*max(times['native16'],times['native32']),reference=2*(fullN+2*pairs*gaussN)/(pilotN+2*gaussN)*times['reference'],comparator=2*fullN/pilotN*times['comparator'],reserve=30)
  targetbytes=sum((data/f'native{N}.toml').stat().st_size for N in c['pilot_meshes']);referencebytes=(data/'reference.json').stat().st_size;rawbytes=targetbytes*fullN/(2*pilotN)+referencebytes*fullN/pilotN;peak=max(rss)+6*rawbytes/1024;forecast=sum(components.values());rc=c['resources'];eligible=max(rss)<=rc['max_pilot_rss_KiB'] and peak<=rc['max_projected_full_rss_KiB'] and forecast<=rc['max_projected_full_seconds']
  for k,v in components.items():close(assessment['components_seconds'][k],v,1e-9,'Independent resource component '+k)
  close(assessment['projected_full_seconds'],forecast,1e-9,'Independent conditional full wall forecast');close(assessment['projected_full_peak_rss_KiB'],peak,1e-6,'Independent conditional full RSS forecast');need(assessment['max_pilot_rss_KiB']==max(rss) and assessment['status']==('ELIGIBLE_FOR_SEPARATE_ROOT_REVIEW' if eligible else 'UNRESOLVED') and assessment['separate_root_current_queue_admission_required'] and not assessment['automatic_submission'],'Exact unchanged engineering decision; no automatic admission')
  report.update(status='PASS',pilot_math_status='PASS',engineering_eligible=eligible,resource_projection=dict(components_seconds=components,full_seconds=forecast,full_peak_RSS_KiB=peak,pilot_max_RSS_KiB=max(rss)),native_rows=8,reference_rows=10,frozen_comparator_gates=len(comparison['checks']),max_recorded_metric_error=max(metrics.values(),default=0.))
 except Exception as ex:report.update(status='UNRESOLVED',reason=type(ex).__name__+': '+str(ex),traceback=traceback.format_exc())
 a.output.mkdir(parents=True);shutil.copyfile(__file__,a.output/'auditor_source.py');inputs={str(p.relative_to(raw)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(raw.rglob('*')) if p.is_file() and p.name!='remote_text_receipts.tar'}
 report.update(schema='J23_INDEPENDENT_EXCLUDED_PILOT_AUDIT_V1',catalogue_id='J23',catalogue_status='NOT_RUN',scientific_completions=0,auditor_sha256=sha(__file__),source_archive_sha256=ARCHIVE,source_freeze_sha256=FREEZE,checks=checks,check_count=len(checks),metric_specific_errors=metrics,inputs=inputs,replay_inputs={p.name:sha(p) for p in sorted(a.replay.iterdir()) if p.is_file()},new_HF_solve_or_quadrature_executed=False,postprocessing_only=True,full_job_admitted=False,scope='Only existing excluded pilot data and prescribed pure trial-function qualification; no SCF/doped-ground/gRISB claim or catalogue completion.')
 (a.output/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');(a.output/'SHA256SUMS').write_text(''.join(sha(p)+'  '+p.name+'\n' for p in sorted(a.output.iterdir()) if p.is_file()));print(json.dumps({k:report[k] for k in ['status','check_count','scientific_completions']}))
if __name__=='__main__':main()
