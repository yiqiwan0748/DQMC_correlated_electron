#!/usr/bin/env python3
"""Independent stdlib-only exact-quarter-character K / sparse Fock-H saved audit."""
import ast,datetime,hashlib,itertools,json,math,pathlib,platform,time,tomllib
from fractions import Fraction
ROOT=pathlib.Path(__file__).resolve().parents[2];SNAP=ROOT/'reconciliation_20261025/weekend_state/snapshot';SRC=SNAP/'source/x19_full_saved_audit_v4_pc1_20261024/auditor';RAW=SNAP/'analysis_inputs/x19_full_packet_6466825_v1';OUT=pathlib.Path(__file__).resolve().parent
start=time.monotonic();checks=[];metrics={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return tomllib.loads(p.read_text())if p.suffix=='.toml'else json.loads(p.read_text())
def need(name,v):
 checks.append(dict(name=name,passed=bool(v)))
 if not v:raise ValueError(name)
def upper(name,v,b):need(name,math.isfinite(v)and 0<=v<=b);metrics[name]=dict(value=v,bound=b)
def defs(p):return{n.name:n for n in ast.parse(p.read_text()).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
def equal(a,b):return ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False)
a=defs(SRC/'audit_x19_full_v4.py');b=defs(SRC/'preserved_v3_source/audit_x19_full_v3.py')
for name in ['sha','load','hop','quotient','exact_torus_K','unpack','independent_H','thermal_observers','replay_equal_only_gate_value','geometry_key','mode_map','open_edges_K','verify_metadata','one_case','map_backends','reconstructed_measure','verify_full_report','qualify_current_maps','mapped_snapshot','loaded_analysis_runtime','actual_imports']:need('unchanged scientific and existing qualification helper '+name,equal(a[name],b[name]))
for name in ['runtime_numpy_scipy_pins.json','current_observed75_registry_PRIVATE.json','optional_cpython_runtime_pins.json']:
 old=SRC/'preserved_v3_source'/name
 if old.exists():need('preserved runtime inventory '+name,sha(SRC/name)==sha(old))
newbody=a['main'].body;oldbody=b['main'].body
# All actual scientific work and frozen comparator replay is inside one Try.
nt=next(v for v in newbody if isinstance(v,ast.Try));ot=next(v for v in oldbody if isinstance(v,ast.Try));need('entire scientific main try/handlers byte AST unchanged',equal(nt,ot))
report=load(SNAP/'results/x19_full_saved_audit_v4_pc1_20261024/audit.json');need('complete full saved audit PASS18462',report['status']=='PASS'and report['full_mathematics']=='PASS'and report['check_count']==len(report['checks'])==18462 and all(v['passed']is True for v in report['checks']));need('full saved matrix audit reached144eigh',report['saved_matrix_math_reached']and report['completed_saved_matrix_eigh_calls']==144);need('strict byte-semantic frozen comparator replay',report['exact_replay_equal']and report['replay_semantics_equal']and report['replay_roundoff_differences']==[])
for name,v in report['metric_specific_errors'].items():upper('saved metric bound '+name,v['value'],v['bound'])
need('all3 guarded failures unchanged',[sha(SRC/'preserved_v3_failure/audit.json'),report['preserved_failed_v2_report_sha256'],report['preserved_failed_v1_report_sha256']]==['be57a92921b4e3d66a74c405a539626372cc901c2c05e714118b581d25534403','918c1a4c74944160c4b6a29ea5e1bca2e37551e16ec36e4c8dee02c988cb816b','ced68791d31ff7a5d362c29760dd71ff3c12bad377bc1ec8edbfe38d295e0a43'])
need('no historical mapping retrofit or producer',not report['historical_failed_runtime_mapping_recovered']and not report['target_reference_producer_or_new_fixture_executed'])
c=load(RAW/'source/manifest.toml');docs=[load(RAW/'results/target.toml'),load(RAW/'results/reference.json')];frozen=load(RAW/'results/comparison.json');replay=load(SNAP/'results/x19_full_saved_audit_v4_pc1_20261024/frozen_comparison_replay.json');need('entire frozen comparison replay exact',frozen==replay)
need('frozen total N4 dimensions',c['N']==4 and c['spin_sectors']==[[0,4],[1,3],[2,2],[3,1],[4,0]]and c['sector_dimensions']==[70,448,784,448,70] and sum(c['sector_dimensions'])==1820)
need('frozen12cases t/beta/h',c['U_arms']==[0.,4.]and c['t4_arms']==[0.,-.65]and(c['beta'],c['h'],c['t1'],c['t2'],c['t3'])==(1.,0.,-1.,-1.75,-.85))
# Exact finite dual characters: all 4cell tori use quarter angles. This uses
# rational roots of unity, not the original floating trigonometric reconstruction.
cos={Fraction(0):1,Fraction(1,4):0,Fraction(1,2):-1,Fraction(3,4):0};sin={Fraction(0):0,Fraction(1,4):1,Fraction(1,2):0,Fraction(3,4):-1}
def trig(q):return cos[q%1],sin[q%1]
def exact_K(case,spec):
 cells=case['cells'];t={k:Fraction(c[k])for k in ['t1','t2','t3']};v=Fraction(case['t4']);K=[[Fraction()for _ in range(8)]for _ in range(8)]
 if spec['kind']=='open':
  need('open exact cellset '+case['id'],{tuple(v)for v in cells}=={(0,0),(0,1),(1,0),(1,1)})
  for i,(x,y)in enumerate(cells):
   for j,(xx,yy)in enumerate(cells):
    dx,dy=xx-x,yy-y
    for orb in [0,1]:
     if abs(dx)==1 and dy==0:K[2*i+orb][2*j+orb]=-t['t1'if orb==0 else't2']
     elif dx==0 and abs(dy)==1:K[2*i+orb][2*j+orb]=-t['t2'if orb==0 else't1']
     elif abs(dx)==abs(dy)==1:K[2*i+orb][2*j+orb]=-t['t3']
    if abs(dx)==abs(dy)==1:K[2*i][2*j+1]=K[2*i+1][2*j]=v*(1 if dx*dy>0 else -1)
 else:
  imaginary=[[Fraction()for _ in range(8)]for _ in range(8)]
  a1,a2=spec['a1'],spec['a2'];D=abs(a1[0]*a2[1]-a1[1]*a2[0]);characters=[(Fraction(p,D),Fraction(q,D))for p in range(D)for q in range(D)if(a1[0]*p+a1[1]*q)%D==0 and(a2[0]*p+a2[1]*q)%D==0];need('all4 exact dual characters '+case['id'],len(characters)==4)
  for px,py in characters:
   cx,sx=trig(px);cy,sy=trig(py);H=[[-2*t['t1']*cx-2*t['t2']*cy-4*t['t3']*cx*cy,-4*v*sx*sy],[-4*v*sx*sy,-2*t['t2']*cx-2*t['t1']*cy-4*t['t3']*cx*cy]]
   for i,(x,y)in enumerate(cells):
    for j,(xx,yy)in enumerate(cells):
     real,imag=trig(px*(x-xx)+py*(y-yy))
     for orb in range(2):
      for other in range(2):
       K[2*i+orb][2*j+other]+=real*H[orb][other]/4;imaginary[2*i+orb][2*j+other]+=imag*H[orb][other]/4
  need('exact Fourier imaginary cancellation '+case['id'],all(v==0 for row in imaginary for v in row))
 return[[float(v)for v in row]for row in K]
def hop(state,i,j):
 if not(state&(1<<j)):return None
 occupied_below_j=(state&((1<<j)-1)).bit_count();reduced=state^(1<<j)
 if reduced&(1<<i):return None
 occupied_below_i=(reduced&((1<<i)-1)).bit_count();return reduced|(1<<i),(-1 if(occupied_below_i+occupied_below_j)%2 else 1)
sector_count=entries=0;maxH=maxK=0.;spectra=[]
for backend,doc in zip(['native','quspin'],docs):
 expected=[(s,U,t4)for s,U,t4 in itertools.product(c['geometries'],c['U_arms'],c['t4_arms'])];need('prescribed12 source cases '+backend,len(doc['cases'])==len(expected)==12)
 for case,(spec,U,t4)in zip(doc['cases'],expected):
  tag=backend+'/'+case['id'];need('prescribed own parameters '+tag,(case['geometry'],case['U'],case['t4'])==(spec['id'],U,t4));K=exact_K(case,spec);err=max(abs(K[i][j]-case['K'][i][j])for i in range(8)for j in range(8));upper('new exact-quarter character K '+tag,err,c['gates']['hamiltonian_abs']);maxK=max(maxK,err);need('all5 prescribed spin sectors '+tag,[[b['Nup'],b['Ndown']]for b in case['sectors']]==c['spin_sectors'])
  energies=[];sector_energies=[]
  for block,dim in zip(case['sectors'],c['sector_dimensions']):
   labels=list(zip(block['basis_up_per_row'],block['basis_down_per_row']));states=[u|(d<<8)for u,d in labels];lookup={s:i for i,s in enumerate(states)};tag2=tag+'/'+str(block['Nup']);need('all physical basis labels complete '+tag2,len(lookup)==dim and set(labels)=={(u,d)for u in range(256)if u.bit_count()==block['Nup']for d in range(256)if d.bit_count()==block['Ndown']});own={}
   for col,(state,(u,d))in enumerate(zip(states,labels)):
    own[col,col]=U*((u&d).bit_count()-(block['Nup']+block['Ndown'])/2+2)
    for spin in [0,1]:
     for i in range(8):
      for j in range(8):
       if K[i][j]:
        moved=hop(state,i+8*spin,j+8*spin)
        if moved:
         dest,sign=moved;key=lookup[dest],col;own[key]=own.get(key,0.)+K[i][j]*sign
   coo=block['H'];saved={(i,j):v for i,j,v in zip(coo['row'],coo['col'],coo['value'])};need('unique finite COO '+tag2,len(saved)==len(coo['row'])==len(coo['col'])==len(coo['value'])and coo['shape']==[dim,dim]and all(0<=i<dim and 0<=j<dim and math.isfinite(v)for(i,j),v in saved.items()));err=max(abs(own.get(k,0.)-saved.get(k,0.))for k in own.keys()|saved.keys());upper('new independent sparse full Fock H '+tag2,err,c['gates']['hamiltonian_abs']);maxH=max(maxH,err);entries+=len(saved);sector_count+=1;need('all saved eigenenergies '+tag2,len(block['energies'])==dim and all(math.isfinite(v)for v in block['energies']));energies.extend(block['energies']);sector_energies.append(block['energies'])
  # Independent canonical weights/logZ/energy directly from all1820 energies.
  e0=min(energies);ws=[math.exp(-c['beta']*(v-e0))for v in energies];Z=math.fsum(ws);logZ=math.log(Z)-c['beta']*e0;energy=math.fsum(v*w for v,w in zip(energies,ws))/Z;o=case['observables'];upper('new canonical logZ '+tag,abs(logZ-o['logZ']),c['gates']['thermal_observable_abs']);upper('new canonical energy '+tag,abs(energy-o['energy']),c['gates']['thermal_observable_abs']);probs=[math.fsum(math.exp(-c['beta']*(e-e0))for e in es)/Z for es in sector_energies];upper('new all5 spin sector probabilities '+tag,max(abs(p-b['canonical_sector_probability'])for p,b in zip(probs,case['sectors'])),c['gates']['normalization_abs']);upper('new total density N4 '+tag,abs(math.fsum(o['density'])-4),c['gates']['normalization_abs']);upper('new canonical total charge covariance '+tag,abs(math.fsum(z for row in o['Cnn']for z in row)),c['gates']['normalization_abs']);upper('new SU2 full correlation '+tag,max(abs(o['Sdot_connected'][i][j]-3*o['Czz'][i][j])for i in range(8)for j in range(8)),c['gates']['spin_rotation_abs'])
need('all120sector sparse matrices reviewed',sector_count==120)
result=dict(schema='WEEKEND_X19_INDEPENDENT_SOURCE_SAVED_DATA_REVIEW_V1',status='PASS',recommendation='Eligible for one X19 catalogue question PASS, subject to root terminal/source/raw/public-packet closure.',scope='Fresh local stdlib-only exact rational quarter-angle dual-character K; complete sparse Fock-H/basis/COO reconstruction for both backends, all120sectors; all1820energies per case canonical weights/energy/logZ and physical sum rules. Scientific v4 helpers and full main scientific try are AST-identical to prior v3. This local analysis does not replay pinned ISAAC runtime or diagonalize new matrices.',python=platform.python_version(),original_pinned_runtime_replay=False,checks=checks,check_count=len(checks),metrics=metrics,sector_blocks=sector_count,raw_COO_entries=entries,max_K_abs=maxK,max_H_abs=maxH,existing_ISAAC_audit_check_count=report['check_count'],existing_ISAAC_audit_sha256=sha(SNAP/'results/x19_full_saved_audit_v4_pc1_20261024/audit.json'),auditor_sha256=sha(SRC/'audit_x19_full_v4.py'),runner_sha256=sha(pathlib.Path(__file__)),inputs={str(p.relative_to(ROOT)):sha(p)for p in[RAW/'source/manifest.toml',RAW/'results/target.toml',RAW/'results/reference.json',RAW/'results/comparison.json',SRC/'audit_x19_full_v4.py',SRC/'preserved_v3_source/audit_x19_full_v3.py']},elapsed_seconds=time.monotonic()-start,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(OUT/'x19_independent_review.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in result.items()if k in ['status','check_count','sector_blocks','raw_COO_entries','max_K_abs','max_H_abs','elapsed_seconds']},indent=2))
