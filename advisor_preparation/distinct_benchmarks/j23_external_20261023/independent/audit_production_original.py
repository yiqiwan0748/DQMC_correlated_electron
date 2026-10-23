#!/usr/bin/env python3
"""Independent saved-data J23 audit: closed projectors, Wick, exact moments, fixed full gates.

No native HF, reference oracle, eigensolver or quadrature is executed. A separate
pinned-runtime frozen-comparator replay is mandatory for a final audit PASS.
"""
import argparse,datetime,hashlib,itertools,json,math,re,shutil,tarfile,tomllib,traceback
from fractions import Fraction
from pathlib import Path
ARCHIVE='5b78ce7ea96a31aedb5ee2c2a27e3f128667ae0121cf410219a13a4ae6ae6bbe'
FREEZE='7b763ad944dd745d604957bf25556fdad06b2ff733f09cb6c25230f7db7b4043'
JOB='6466733'
PILOT_AUDITOR='899360edd906a42ce6293a6b595dd9788314cdc41bbd6929907c6659b02b7416'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return tomllib.loads(p.read_text()) if p.suffix=='.toml' else json.loads(p.read_text())
def flat(x):
 if isinstance(x,list):
  for v in x:yield from flat(v)
 else:yield x
def shape(x):
 if not isinstance(x,list):return ()
 sizes=[shape(v) for v in x]
 if sizes and any(s!=sizes[0] for s in sizes):raise ValueError('Ragged scientific array')
 return (len(x),)+(sizes[0] if sizes else ())
def maxerr(a,b):
 if shape(a)!=shape(b):raise ValueError('Observable shape mismatch')
 aa=list(flat(a));bb=list(flat(b))
 if not all(type(x) in (float,int) and math.isfinite(x) for x in aa+bb):raise ValueError('Nonfinite/nonreal observable')
 return max((abs(x-y) for x,y in zip(aa,bb)),default=0.)
def sub(a,b):
 if shape(a)!=shape(b):raise ValueError('Difference shape mismatch')
 return [sub(x,y) for x,y in zip(a,b)] if isinstance(a,list) else a-b

def independent_native(row,c):
 """Own closed2x2 projectors on saved N, including exactly degenerate filled blocks."""
 N=row['mesh'];m=row['m'];t4=row['t4'];V=[[-m,m],[-m,m]] if row['alignment']=='parallel' else [[-m,m],[m,-m]]
 rho=[[0.,0.],[0.,0.]];kin=0.;bloch=[];eigs=[[],[]];occ=[[],[]];gap=math.inf;counts=[]
 for i in range(N):
  x=-math.pi+2*math.pi*(i+c['native_mesh_shift'])/N
  for j in range(N):
   y=-math.pi+2*math.pi*(j+c['native_mesh_shift'])/N
   ex=-2*c['t1']*math.cos(x)-2*c['t2']*math.cos(y)-4*c['t3']*math.cos(x)*math.cos(y)
   ey=-2*c['t2']*math.cos(x)-2*c['t1']*math.cos(y)-4*c['t3']*math.cos(x)*math.cos(y)
   h=-4*t4*math.sin(x)*math.sin(y);bloch.append([ex,ey,h,0.]);count=0
   for spin in range(2):
    a=(ex+V[0][spin]+ey+V[1][spin])/2;z=(ex+V[0][spin]-ey-V[1][spin])/2;R=math.hypot(z,h)
    energies=[a-R,a+R];f=[int(e<row['native_evaluate']['mu']) for e in energies];eigs[spin].append(energies);occ[spin].append(f);count+=sum(f);gap=min(gap,*(abs(e-row['native_evaluate']['mu']) for e in energies))
    if R==0:
     if f[0]!=f[1]:raise ValueError('Unqualified partial occupation of degenerate2x2 block')
     rho[0][spin]+=f[0]/(N*N);rho[1][spin]+=f[0]/(N*N);kin+=f[0]*(ex+ey)/(N*N)
    else:
     for sign,occupied in zip([-1,1],f):
      if occupied:
       p0=(1+sign*z/R)/2;p1=(1-sign*z/R)/2;p01=sign*h/(2*R)
       rho[0][spin]+=p0/(N*N);rho[1][spin]+=p1/(N*N);kin+=(ex*p0+ey*p1+2*h*p01)/(N*N)
   counts.append(count)
 F=kin+2*m*sum((r[0]-.5)*(r[1]-.5) for r in rho)
 return dict(rho=rho,kinetic=kin,physicalF=F,V=V,bloch=bloch,eigenvalues=eigs,occupations=occ,gap=gap,counts=counts)

def exact_coefficient(c,m,nmax=8):
 """No quadrature: exact cosine moments and a geometric-series tail enclosure."""
 def moment(k):return Fraction(0) if k%2 else Fraction(math.comb(k,k//2),2**k)
 def A(k):return moment(k)-moment(k+2)
 mf=Fraction(str(m));af=Fraction(str(c['t2']))-Fraction(str(c['t1']));ratio=af*af/(mf*mf);q=4*ratio
 if not 0<=q<1:raise ValueError('Analytic moment series requires |2(t2-t1)/m|<1')
 s=Fraction(0)
 for n in range(nmax+1):
  Mn=sum((Fraction(math.comb(2*n,j))*A(2*n-j)*A(j) for j in range(0,2*n+1,2)),Fraction(0))
  s+=ratio**n*Mn
 estimate=-16*s/mf;tail=4*q**(nmax+1)/(mf*(1-q))
 # Float enclosure is outward rounded and includes the conversion of both rationals.
 return dict(center=float(estimate),tail=math.nextafter(float(tail)+8*math.ulp(float(estimate)),math.inf),q=str(q),terms=nmax+1,exact_center=str(estimate),exact_tail=str(tail))

class FixedGates:
 """Independent stdlib reconstruction of the original registered full gates."""
 def __init__(self):self.rows=[]
 def check(self,name,value,relation,bound,kind='qualification'):
  if relation=='eq':passed=value==bound
  else:passed=math.isfinite(value) and (value<=bound if relation=='le' else value>=bound if relation=='ge' else value<bound)
  self.rows.append(dict(name=name,value=value,relation=relation,bound=bound,kind=kind,passed=bool(passed)))
 def error(self,name,x,y,limit,kind='comparison'):self.check(name,maxerr(x,y),'le',limit,kind)

def rebuild_fixed_gates(native,reference,c):
 g=FixedGates();L=c['gates'];params=list(itertools.product(c['m_values'],c['t4_values'],c['alignments']))
 refs={kind:{(r['size'],r['m'],r['t4'],r['alignment']):r for r in reference['rows'] if r['grid_kind']==kind} for kind in ['native_mesh','gauss']}
 for kind,sizes in [('native_mesh',c['native_square_mesh']),('gauss',c['oracle_gauss_orders'])]:
  exp={(N,m,t4,a) for N in sizes for m,t4,a in params}
  if set(refs[kind])!=exp or len(refs[kind])!=sum(r['grid_kind']==kind for r in reference['rows']):raise ValueError('Missing/duplicate full reference grid')
 table={(r['mesh'],r['m'],r['t4'],r['alignment'],r['phase']):r for r in native['rows']};expected={(N,m,t4,a,'production') for N in c['native_square_mesh'] for m,t4,a in params}
 if len(native['rows'])!=len(table) or set(table)!=expected:raise ValueError('Missing/duplicate full native grid')
 for (N,m,t4,a,phase),r in table.items():
  tag=f'N{N}/m{m}/t4{t4}/{a}/{phase}';ref=refs['native_mesh'][N,m,t4,a];actual=r['native_evaluate'];obs=r['observer'];U=2*m;rho0=[[1.,0.],[1.,0.]] if a=='parallel' else [[1.,0.],[0.,1.]]
  if r['Ns']!=1 or r['d']!=2 or r['nk']!=N*N or r['U']!=U:raise ValueError('Native normalization metadata')
  g.error(tag+'/actual trial input',r['rho0'],rho0,L['direct_analytic_eigenvalues_abs'],'qualification');g.error(tag+'/actual native V',actual['V'],[[U*(z[1]-.5),U*(z[0]-.5)] for z in rho0],L['direct_analytic_eigenvalues_abs'],'qualification')
  eig=obs['eigenvalues_by_spin_k_band'];f=obs['occupations_by_spin_k_band']
  if shape(eig)!=(2,N*N,2) or shape(f)!=shape(eig) or not all(math.isfinite(v) for v in itertools.chain(flat(eig),flat(f))):raise ValueError('Complete eigen/occupation shape')
  counts=[sum(f[s][k][b] for s in range(2) for b in range(2)) for k in range(N*N)]
  g.check(tag+'/binary pure occupations',all(v in (0,1) for v in flat(f)) and obs['pure_exact_binary_occupations'],'eq',True);g.check(tag+'/two occupied at every k',all(v==2 for v in counts) and obs['two_occupied_per_k'],'eq',True);g.error(tag+'/saved occupied count',obs['occupations_per_k'],[2.]*(N*N),0.,'qualification')
  occ=[actual['mu']-en for en,ff in zip(flat(eig),flat(f)) if ff==1];empty=[en-actual['mu'] for en,ff in zip(flat(eig),flat(f)) if ff==0];gap=min(min(occ),min(empty))
  g.check(tag+'/actual common mu gap',gap,'ge',L['common_mu_gap_min']);g.check(tag+'/saved common mu gap',abs(gap-obs['common_mu_gap_min']),'le',L['direct_analytic_eigenvalues_abs']);g.check(tag+'/actual density',abs(actual['n']-2),'le',L['native_density_abs']);g.check(tag+'/actual entropy',abs(actual['entropy']),'le',L['native_entropy_abs']);g.check(tag+'/observer entropy',abs(obs['entropy']),'le',L['native_entropy_abs']);g.check(tag+'/eigenpair residual',obs['maximum_eigenpair_residual'],'le',L['direct_analytic_eigenvalues_abs'])
  g.error(tag+'/observer reproduces actual rho',obs['rho'],actual['rho'],L['observer_density_vs_actual_native_rho_abs'],'qualification');g.check(tag+'/observer reproduces actual physicalF',abs(obs['physicalF']-actual['physicalF']),'le',L['same_mesh_physical_energy_per_cell_abs']);g.check(tag+'/observer reproduces actual kinetic',abs(obs['kinetic']-actual['kinetic']),'le',L['same_mesh_physical_energy_per_cell_abs'])
  g.error(tag+'/all actual Bloch entries',obs['bloch_entries_by_k'],ref['bloch_entries_by_k'],L['direct_analytic_eigenvalues_abs']);g.error(tag+'/all analytic AppendixD eigenvalues',eig,ref['eigenvalues_by_spin_k_band'],L['direct_analytic_eigenvalues_abs']);g.error(tag+'/independent pure occupation prescription',f,ref['occupations_by_spin_k_band'],0.);g.error(tag+'/same mesh actual rho',actual['rho'],ref['rho'],L['observer_density_vs_actual_native_rho_abs'])
  for field in ['physicalF','kinetic']:g.check(tag+'/same mesh independent physical '+field,abs(actual[field]-ref[field]),'le',L['same_mesh_physical_energy_per_cell_abs'],'comparison')
 for kind,lookup in refs.items():
  for k,r in lookup.items():
   tag=kind+'/'+str(k);g.check(tag+'/analytic pure-state qualification',r['pure_occupations_qualified'],'eq',True);g.check(tag+'/common gap',r['common_mu_gap_min'],'ge',L['common_mu_gap_min']);g.check(tag+'/density',abs(r['density']-2),'le',L['native_density_abs']);g.check(tag+'/integration measure normalization',abs(r['quadrature_weight_sum']-1),'le',L['native_density_abs'])
 nat={(N,m,t4,a):r['native_evaluate'] for (N,m,t4,a,_),r in table.items()}
 def difference(lookup,N,m,t4,field):return sub(lookup[N,m,t4,'antiparallel'][field],lookup[N,m,t4,'parallel'][field])
 for m,t4 in itertools.product(c['m_values'],c['t4_values']):
  for field in ['physicalF','kinetic','rho']:
   for a in c['alignments']:
    g.error(f'native64to128/m{m}/t4{t4}/{a}/{field}',nat[64,m,t4,a][field],nat[128,m,t4,a][field],L['final_native_mesh_difference_abs'],'qualification');g.error(f'Gauss96to128/m{m}/t4{t4}/{a}/{field}',refs['gauss'][96,m,t4,a][field],refs['gauss'][128,m,t4,a][field],L['oracle_quadrature_difference_abs'],'qualification');g.error(f'final independent/m{m}/t4{t4}/{a}/{field}',nat[128,m,t4,a][field],refs['gauss'][128,m,t4,a][field],L['observer_density_vs_actual_native_rho_abs'] if field=='rho' else L['native128_vs_oracle128_energy_abs'])
   g.error(f'native delta stability/m{m}/t4{t4}/{field}',difference(nat,64,m,t4,field),difference(nat,128,m,t4,field),L['final_native_mesh_difference_abs'],'qualification');g.error(f'Gauss delta stability/m{m}/t4{t4}/{field}',difference(refs['gauss'],96,m,t4,field),difference(refs['gauss'],128,m,t4,field),L['oracle_quadrature_difference_abs'],'qualification');g.error(f'final delta independent/m{m}/t4{t4}/{field}',difference(nat,128,m,t4,field),difference(refs['gauss'],128,m,t4,field),L['observer_density_vs_actual_native_rho_abs'] if field=='rho' else L['native128_vs_oracle128_energy_abs'])
 for name,lookup,sizes in [('native',nat,c['native_square_mesh']),('Gauss',refs['gauss'],c['oracle_gauss_orders'])]:
  for N,m in itertools.product(sizes,c['m_values']):
   g.check(f'{name}{N}/m{m}/zero hybrid alignment',abs(difference(lookup,N,m,0.,'physicalF')),'le',L['zero_t4_alignment_energy_difference_abs'])
   for a in c['alignments']:g.check(f'{name}{N}/m{m}/{a}/even t4',abs(lookup[N,m,.04,a]['physicalF']-lookup[N,m,-.04,a]['physicalF']),'le',L['even_t4_energy_abs'])
   for t4 in c['t4_values']:
    if t4:g.check(f'{name}{N}/m{m}/t4{t4}/antiparallel lower physical trial energy',difference(lookup,N,m,t4,'physicalF'),'le',L['nonzero_hybrid_delta_energy_upper'])
 details=[]
 for m in c['m_values']:
  analytic=refs['gauss'][128,m,0.,'parallel']['leading_coefficient'];g.check(f'm{m}/leading coefficient Gauss96to128',abs(analytic-refs['gauss'][96,m,0.,'parallel']['leading_coefficient']),'le',L['oracle_quadrature_difference_abs']);estimates={}
  for backend,lookup in [('native',nat),('Gauss',refs['gauss'])]:
   zero=difference(lookup,128,m,0.,'physicalF')
   def Q(h):return (difference(lookup,128,m,h,'physicalF')-zero)/(h*h)
   R1=(4*Q(.01)-Q(.02))/3;R2=(4*Q(.02)-Q(.04))/3;estimates[backend]=dict(primary=R1,heldout=R2)
   g.check(f'{backend}/m{m}/primary Richardson vs analytic coefficient',abs(R1-analytic),'le',L['leading_coefficient_abs'],'comparison');g.check(f'{backend}/m{m}/heldout Richardson vs analytic coefficient',abs(R2-analytic),'le',L['leading_coefficient_abs'],'comparison');g.check(f'{backend}/m{m}/primary vs heldout no selection',abs(R1-R2),'le',L['Richardson_primary_vs_heldout_abs'],'qualification')
  details.append(dict(m=m,analytic_leading_coefficient=analytic,Richardson=estimates))
 quality=any(not x['passed'] and x['kind']=='qualification' for x in g.rows);bad=any(not x['passed'] for x in g.rows)
 return g.rows,details,'UNRESOLVED' if quality else 'FAIL' if bad else 'PASS'

def equivalent(x,y,path=(),changes=None):
 """Only declared numeric gate-value roundoff may vary; decisions/bounds exact."""
 if type(x) is not type(y):return False
 if isinstance(x,dict):return x.keys()==y.keys() and all(equivalent(v,y[k],path+(k,),changes) for k,v in x.items())
 if isinstance(x,list):return len(x)==len(y) and all(equivalent(v,w,path+(i,),changes) for i,(v,w) in enumerate(zip(x,y)))
 if x==y:return True
 if isinstance(x,float) and len(path)>=3 and path[0]=='checks' and isinstance(path[1],int) and path[-1]=='value' and math.isfinite(x) and math.isfinite(y):
  tol=max(1e-12,64*math.ulp(1.)*max(1.,abs(x),abs(y)))
  if abs(x-y)<=tol:
   if changes is not None:changes.append(dict(path=list(path),original=x,replayed=y,absolute_difference=abs(x-y),tolerance=tol))
   return True
 return False

def main():
 p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--replay',type=Path,required=True,help='New compute-allocation replay folder containing comparison_replay.json');p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists();source=a.raw/'source';data=a.raw/'results';checks=[];metrics={};report={};physics=[];replay_changes=[]
 def need(v,n):checks.append(dict(name=n,passed=bool(v)));assert v,n
 def close(x,y,limit,n,kind='comparison'):
  value=maxerr(x,y);ok=value<=limit;metrics[n]=value;physics.append(dict(name=n,value=value,bound=limit,kind=kind,passed=ok))
 try:
  fr=load(a.raw/'release_freeze.json');c=load(source/'manifest.toml');ex=load(data/'execution.json');comparison=load(data/'comparison.json');reference=load(data/'reference.json');native=load(data/'native.toml');fresh=load(a.replay/'comparison_replay.json')
  need(sha(a.raw/'release_freeze.json')==FREEZE and sha(a.raw/'source_archive.tar.gz')==fr['archive_sha256']==ARCHIVE,'Exact original frozen full source/archive')
  with tarfile.open(a.raw/'source_archive.tar.gz') as tf:
   ms=tf.getmembers();need(len(ms)==len(fr['files'])==19,'All19 exact frozen source members')
   for m in ms:
    need(m.isfile() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts,'Safe regular archive member');n=str(Path(m.name).relative_to(fr['prefix'].rstrip('/')));b=tf.extractfile(m).read();need(hashlib.sha256(b).hexdigest()==fr['files'][n] and b==(source/n).read_bytes(),'Exact frozen source '+n)
  need(sha(source/'source_freeze.json')==FREEZE,'Source-local original freeze');sourcehash={n:h for h,n in (line.split('  ',1) for line in (source/'SHA256SUMS').read_text().splitlines())};need(len(sourcehash)==18 and all(sha(source/n)==h for n,h in sourcehash.items()),'Exact18 substantive source hashes')
  before=load(a.raw/'transfer_before.json');after=load(a.raw/'transfer_after.json');amap={r['local']:r for r in after['files']};need(before['job_id']==JOB and after['status']=='PASS','Actual full-job unchanged remote-after receipt')
  for r in before['files']:
   q=amap[r['local']];f=a.raw/r['local'];need(f.is_file() and not f.is_symlink() and f.stat().st_size==r['bytes']==q['bytes'] and sha(f)==r['sha256']==q['sha256'] and q['mtime_ns']==r['before']['mtime_ns'] and q['inode']==r['before']['inode'],'Transferred original hash/bytes/stat '+r['local'])
  lib=load(source/'reference_library_pins.json')
  for n,h in lib['files'].items():
   r=before['runtime_files']['reference/'+n];q=after['runtime_files']['reference/'+n];need(r['sha256']==q['sha256']==h and r['bytes']==q['bytes']==lib['bytes'][n] and r['mtime_ns']==q['mtime_ns'] and r['inode']==q['inode'],'Actual39-file source/native/reference pin '+n)
  for role in ['Python_executable','Julia_executable']:
   r=before['runtime_files'][role];q=after['runtime_files'][role];need(r['sha256']==q['sha256'] and r['bytes']==q['bytes'] and r['mtime_ns']==q['mtime_ns'] and r['inode']==q['inode'],'Actual unchanged executable '+role)
  jobs=[s.split('|') for s in (a.raw/'sacct_fresh.txt').read_text().splitlines() if s.startswith(JOB+'|')];need(len(jobs)==1 and jobs[0][2:4]==['COMPLETED','0:0'] and jobs[0][5:8]==['1','3500M','00:30:00'],'Original full accounting1CPU3500MB30minutes')
  need(ex['schema']=='J23_EXECUTION_RECEIPT_V1' and ex['mode']=='production' and ex['job_id']==JOB and ex['source_freeze_sha256']==FREEZE and ex['source_sha256']==sourcehash and ex['scientific_completions']==0 and not ex['SCF_executed'],'Actual source-bound original full execution')
  need(ex['cpus']==1 and ex['memory_MB']==3500 and ex['finished_unix']-ex['started_unix']<=c['resources']['production_combined_timeout_seconds'],'Fixed full CPU/memory/1620s bound')
  need([x['name'] for x in ex['processes']]==['native','reference','comparator'] and all(not x['timed_out'] and x['elapsed_seconds']>0 for x in ex['processes']),'Complete3process full grid, no timeout')
  # Frozen run_campaign treats comparator FAIL/U as an incomplete process; retain that original status rather than relabel it.
  need(all(x['exit_code']==0 for x in ex['processes'][:2]) and ex['processes'][2]['exit_code'] in (0,1,2),'Both actual producers complete; original comparator exit retained')
  for n,r in ex['outputs'].items():f=data/n;need(f.is_file() and sha(f)==r['sha256'] and f.stat().st_size==r['bytes'],'Execution-bound output '+n)
  need(native['schema']=='J23_ACTUAL_NATIVE_HF_RAW_V1' and native['mode']=='production' and native['job_id']==JOB and native['source_sha256']==sourcehash and native['source_freeze_sha256']==FREEZE and native['native_source_sha256']==c['native_source_hashes'] and native['julia_version']==c['julia_version'] and native['trial_state_only'] and not native['SCF_executed'] and native['Ns']==1 and native['d']==2,'Actual unchanged native APIs/version/normalization/trial-only scope')
  need(reference['schema']=='J23_INDEPENDENT_ANALYTIC_WICK_RAW_V1' and reference['mode']=='production' and reference['job_id']==JOB and reference['source_sha256']==sourcehash and reference['source_freeze_sha256']==FREEZE and reference['reference_library_pins_sha256']==sha(source/'reference_library_pins.json') and reference['reference_library_sha256']==lib['files'],'Exact independent reference closure')
  need(reference['python_version']==c['python_version'] and reference['versions']=={'numpy':c['numpy_version'],'scipy':c['scipy_version']} and not reference['SCF_executed'] and not reference['ground_state_claim'],'Pinned actual reference runtime and limited trial scope')
  expected={(N,m,t,a,'production') for N,m,t,a in itertools.product(c['native_square_mesh'],c['m_values'],c['t4_values'],c['alignments'])};actual={(r['mesh'],r['m'],r['t4'],r['alignment'],r['phase']) for r in native['rows']};need(len(native['rows'])==len(actual)==112 and actual==expected,'All112 distinct full native rows')
  grids={(r['mesh'],r['t4']) for r in native['grids']};need(len(native['grids'])==len(grids)==28 and grids==set(itertools.product(c['native_square_mesh'],c['t4_values'])),'All28 actual native Bloch grids')
  refs={(r['grid_kind'],r['size'],r['m'],r['t4'],r['alignment']):r for r in reference['rows']};rexpect={(kind,N,m,t,a) for kind,ns in [('native_mesh',c['native_square_mesh']),('gauss',c['oracle_gauss_orders'])] for N,m,t,a in itertools.product(ns,c['m_values'],c['t4_values'],c['alignments'])};need(len(reference['rows'])==len(refs)==196 and set(refs)==rexpect,'All196 reference rows:112native-mesh and84Gauss')
  coefficients={m:exact_coefficient(c,m) for m in c['m_values']}
  for r in native['rows']:
   tag=f"N{r['mesh']}/m{r['m']}/t4{r['t4']}/{r['alignment']}";own=independent_native(r,c);obs=r['observer'];actual=r['native_evaluate'];L=c['gates'];ref=refs['native_mesh',r['mesh'],r['m'],r['t4'],r['alignment']]
   close(obs['bloch_entries_by_k'],own['bloch'],L['direct_analytic_eigenvalues_abs'],'Independent complete Bloch '+tag);close(obs['eigenvalues_by_spin_k_band'],own['eigenvalues'],L['direct_analytic_eigenvalues_abs'],'Independent complete eigenvalues '+tag);close(obs['occupations_by_spin_k_band'],own['occupations'],0.,'Independent exact occupations '+tag,'qualification');close(actual['V'],own['V'],L['direct_analytic_eigenvalues_abs'],'Independent physical exchange signs '+tag,'qualification');close(obs['occupations_per_k'],own['counts'],0.,'Independent every-k count '+tag,'qualification');need(own['gap']>=L['common_mu_gap_min'] and all(z==2 for z in own['counts']),'Own closed-projector gapped pure qualification '+tag)
   for field in ['rho','kinetic','physicalF']:
    tol=L['observer_density_vs_actual_native_rho_abs'] if field=='rho' else L['same_mesh_physical_energy_per_cell_abs'];close(actual[field],own[field],tol,'Independent physical actual '+field+'/'+tag);close(ref[field],own[field],tol,'Independent saved reference '+field+'/'+tag)
   close(actual['n'],2.,L['native_density_abs'],'Independent density '+tag,'qualification');close(actual['entropy'],0.,L['native_entropy_abs'],'Independent zero entropy '+tag,'qualification');close(obs['rho'],actual['rho'],L['observer_density_vs_actual_native_rho_abs'],'Observer qualifies against actual rho '+tag,'qualification')
  for key,r in refs.items():
   tag=str(key);need(r['nk']==r['size']**2 and r['U']==2*r['m'] and r['pure_occupations_qualified'] and r['common_mu_gap_min']>=1.,'Reference fixed pure/gapped physical metadata '+tag);close(r['quadrature_weight_sum'],1.,1e-12,'Reference measure '+tag,'qualification');close(r['density'],2.,1e-12,'Reference density '+tag,'qualification');close(r['entropy'],0.,1e-12,'Reference entropy '+tag,'qualification');need(all(0<=v<=1 for v in flat(r['rho'])),'Reference positive local occupation '+tag)
   wick=r['kinetic']+r['U']*sum((p[0]-.5)*(p[1]-.5) for p in r['rho']);close(r['physicalF'],wick,1e-10,'Independent reference Wick '+tag)
   cert=coefficients[r['m']];close(r['leading_coefficient'],cert['center'],c['gates']['oracle_quadrature_difference_abs']+cert['tail'],'Independent exact-cosine-moment coefficient '+tag)
   if r['alignment']=='parallel':close(r['rho'],[[1.,0.],[1.,0.]],1e-12,'Parallel filled spin projector '+tag);close(r['kinetic'],0.,1e-10,'Parallel complete-band kinetic '+tag);close(r['physicalF'],-r['m'],1e-10,'Parallel analytic physical functional '+tag)
  gates,details,status=rebuild_fixed_gates(native,reference,c);need(len(gates)==len(comparison['checks']),'Every frozen gate independently reconstructed')
  need(equivalent({'checks':comparison['checks']},{'checks':gates},changes=replay_changes),'Independent reconstruction: every gate name,bound,kind,decision and value')
  need(comparison['native_rows']==112 and comparison['reference_rows']==196 and comparison['mode']=='production' and not comparison['pilot_excluded'] and comparison['source_freeze_sha256']==FREEZE and comparison['input_sha256']=={n:sha(data/n) for n in ['native.toml','reference.json']},'Original complete comparator grid/input bindings')
  need(comparison['status']==comparison['catalogue_status']==status and comparison['scientific_completions']==int(status=='PASS'),'Original science classification agrees with all fixed gates')
  need(maxerr([d['analytic_leading_coefficient'] for d in details],[d['analytic_leading_coefficient'] for d in comparison['leading_coefficients']])==0.,'Exact independently rebuilt coefficient anchors')
  need(equivalent(comparison,fresh,changes=replay_changes),'Fresh unchanged pinned-runtime frozen comparator replay; only gate-value arithmetic tolerance')
  ownbad=[x for x in physics if not x['passed']];need(not ownbad,'Independent closed-projector/Wick/exact-moment checks agree')
  report.update(status=status,audit_integrity='PASS',science_status=status,catalogue_status=status,scientific_completions=int(status=='PASS'),native_rows=112,reference_rows=196,Gauss_reference_rows=84,independently_rebuilt_frozen_gates=len(gates),frozen_gate_failures=[r for r in gates if not r['passed']],leading_coefficients=details,independent_exact_moment_certificates=coefficients,original_execution_status=ex['status'])
 except Exception as err:report.update(status='UNRESOLVED',audit_integrity='UNRESOLVED',science_status='UNRESOLVED',catalogue_status='UNRESOLVED',scientific_completions=0,reason=type(err).__name__+': '+str(err),traceback=traceback.format_exc())
 a.output.mkdir(parents=True);shutil.copyfile(__file__,a.output/'auditor_source.py');inputs={str(p.relative_to(a.raw)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.raw.rglob('*')) if p.is_file() and p.name!='remote_text_receipts.tar'}
 report.update(schema='J23_INDEPENDENT_FULL_SAVED_DATA_AUDIT_V1',catalogue_ids=['J23'],job_id=JOB,auditor_sha256=sha(__file__),source_archive_sha256=ARCHIVE,source_freeze_sha256=FREEZE,pilot_mathematical_predecessor_sha256=PILOT_AUDITOR,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=checks,check_count=len(checks),independent_physics_checks=physics,metric_specific_errors=metrics,portable_gate_value_differences=replay_changes,inputs=inputs,input_sha256={k:v['sha256'] for k,v in inputs.items()},replay_inputs={str(p.relative_to(a.replay)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in a.replay.rglob('*') if p.is_file()},new_HF_or_quadrature_or_eigensolver_executed=False,postprocessing_only=True,scope='One prescribed n2 pure-Slater trial functional at actual returned native fields; no SCF,doped phase,gRISB or exact-ground-state claim. Independent closed-form projectors and exact cosine moments use saved model/grids; no new scientific producer.')
 (a.output/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');(a.output/'SHA256SUMS').write_text(''.join(sha(p)+'  '+p.name+'\n' for p in sorted(a.output.iterdir()) if p.is_file()));print(json.dumps({k:report[k] for k in ['status','audit_integrity','check_count','scientific_completions']}))
if __name__=='__main__':main()
