#!/usr/bin/env python3
"""Audit retained eigensystems/operators via independent Fock algebra; no solver rerun."""
import argparse,hashlib,itertools,json,math,tomllib
from pathlib import Path
import numpy as np
FREEZE='469eda04dfaf88287513fdf670871057cb0f4c24a65fd996ce8be9268ba58bc2'
ARCHIVE='a92f18307173a500317563ed888ccc5c3d1b44c53ad8d6616a77a52e3558074f'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return tomllib.loads(p.read_text()) if p.suffix=='.toml' else json.loads(p.read_text())
def decode(x):
 a=np.zeros(x['shape']);i,j=x['row'],x['col'];assert len(set(zip(i,j)))==len(i)==len(x['value']);a[i,j]=x['value'];return a

def hop(f,p,q):
 if not(f>>q&1):return None
 s=-1 if (f&((1<<q)-1)).bit_count()%2 else 1;f^=1<<q
 if f>>p&1:return None
 s*= -1 if (f&((1<<p)-1)).bit_count()%2 else 1
 return f|(1<<p),s

def geom(cfg,spec):
 nx,ny=spec.get('Lx',cfg['Lx']),spec.get('Ly',cfg['Ly']);M=2*nx*ny;K=np.zeros((M,M))
 # Independent coordinate-pair geometry, rather than production bond loop/Kronecker code.
 modes=[(x,y,o) for y in range(ny) for x in range(nx) for o in range(2)]
 for i,(x,y,o) in enumerate(modes):
  for j,(a,b,c) in enumerate(modes):
   dx,dy=a-x,b-y
   if o==c and abs(dx)==1 and dy==0:K[i,j]=-cfg['t1' if o==0 else 't2']
   if o==c and dx==0 and abs(dy)==1:K[i,j]=-cfg['t2' if o==0 else 't1']
   if abs(dx)==abs(dy)==1:K[i,j]=-cfg['t3'] if o==c else cfg['t4']*dx*dy
 return K

def independent_operators(coords,K,U):
 M=len(K);dim=len(coords);focks=[u|(d<<M) for u,d in coords];idx={f:i for i,f in enumerate(focks)};H=np.zeros((dim,dim));C=np.zeros_like(H)
 for j,(u,d) in enumerate(coords):
  H[j,j]=U*(u&d).bit_count()
  for sigma in (0,1):
   for p,q in zip(*np.nonzero(K)):
    z=hop(focks[j],int(p)+sigma*M,int(q)+sigma*M)
    if z:H[idx[z[0]],j]+=K[p,q]*z[1]
   for p in range(0,M,2):
    z=hop(focks[j],p+sigma*M,p+1+sigma*M)
    if z:C[idx[z[0]],j]+=z[1]
 nu,nd=coords[0][0].bit_count(),coords[0][1].bit_count();Sp={}
 for j,f in enumerate(focks):
  for k in range(M):
   z=hop(f,k,k+M)
   if z:
    if z[0] not in Sp:Sp[z[0]]=np.zeros(dim)
    Sp[z[0]][j]+=z[1]
 B=np.array(list(Sp.values()));S=B.T@B if len(B) else np.zeros((dim,dim));sz=(nu-nd)/2;S+=sz*(sz+1)*np.eye(dim)
 return H,C,S

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--receipt-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--public-input-manifest',type=Path);ap.add_argument('--public-execution-metadata',type=Path);ap.add_argument('--expected-public-manifest-sha256');ap.add_argument('--expected-original-execution-sha256');a=ap.parse_args();p=a.receipt_root;assert not a.output.exists();checks=[];metrics=[]
 def need(n,b):checks.append({'name':n,'passed':bool(b)});assert b,n
 def close(n,x,y,tol):
  e=float(np.max(np.abs(np.asarray(x)-np.asarray(y)))) if np.size(x) else 0.;checks.append({'name':n,'error':e,'limit':tol,'passed':math.isfinite(e) and e<=tol});assert checks[-1]['passed'],n;return e
 source=p/'source';cfg=load(source/'manifest.toml');freeze=load(source/'source_freeze.json');need('exact registered source freeze',sha(source/'source_freeze.json')==FREEZE and freeze['archive_sha256']==ARCHIVE)
 for n,h in freeze['files'].items():need('source '+n,sha(source/n)==h)
 for line in (source/'SHA256SUMS').read_text().splitlines():
  h,n=line.split('  ',1);need('producer source pin '+n,sha(source/n)==h)
 public=bool(a.public_input_manifest or a.public_execution_metadata)
 if public:
  need('public mode explicit integrity arguments',bool(a.public_input_manifest and a.expected_public_manifest_sha256));need('pinned public input manifest',sha(a.public_input_manifest)==a.expected_public_manifest_sha256);pm=load(a.public_input_manifest);need('public manifest schema',pm['schema']=='PUBLIC_PYSCF_SCIENTIFIC_INPUT_MANIFEST_V1')
  for n,r in pm['files'].items():need('public input hash '+n,sha(p/n)==r['sha256'] and (p/n).stat().st_size==r['bytes'])
  need('public manifest covers scientific raw and source freeze',all(n in pm['files'] for n in ['results/target.toml','results/reference.json','results/comparison.json','source/source_freeze.json']))
  if a.public_execution_metadata:
   execution=load(a.public_execution_metadata);need('explicit public execution projection',bool(a.expected_original_execution_sha256) and execution['schema']=='PUBLIC_PYSCF_EXECUTION_PROJECTION_V1' and execution['original_schema']=='PYSCF_TWO_ORBITAL_EXECUTION_V1' and execution['original_execution_sha256']==a.expected_original_execution_sha256);execution_sha=a.expected_original_execution_sha256
  else:
   need('public original execution covered', 'results/execution.json' in pm['files']);execution=load(p/'results/execution.json');execution_sha=sha(p/'results/execution.json')
 else:
  before=load(p/'transfer_before.json');after=load(p/'transfer_after.json');need('remote after verified',after['status']=='PASS');aft={r['local']:r for r in after['files']}
  for r in before['files']:
   need('input hash '+r['local'],sha(p/r['local'])==r['sha256'] and (p/r['local']).stat().st_size==r['bytes']);z=aft[r['local']];need('remote unchanged '+r['local'],z['sha256']==r['sha256'] and z['bytes']==r['bytes'] and z['mtime_ns']==r['before']['mtime_ns'] and z['inode']==r['before']['inode'])
  execution=load(p/'results/execution.json');execution_sha=sha(p/'results/execution.json')
 raw={k:load(p/'results'/n) for k,n in [('target','target.toml'),('reference','reference.json')]};comp=load(p/'results/comparison.json');replay=load(p/'independent_replay/comparison.json')
 need('production all processes complete',execution['status']=='PROCESSES_COMPLETE' and execution['mode']=='production' and execution['cpu_limit']==1 and all(r['exit_code']==0 for r in execution['processes']) and execution['source_freeze_sha256']==FREEZE)
 need('workflow cap',execution['finished_unix']<=execution['deadline_unix'] and execution['finished_unix']-execution['started_unix']<=cfg['production_timeout_seconds'])
 for n,r in execution['outputs'].items():need('execution raw hash '+n,sha(p/'results'/n)==r['sha256'] and (p/'results'/n).stat().st_size==r['bytes'])
 expected=[s['id'] for s in cfg['cases']];need('all7 registered fixtures',len(expected)==7 and all([r['id'] for r in x['cases']]==expected and x['mode']=='production' and x['resource_repeat_cases']==[] and x['threads']==1 for x in raw.values()))
 pins=load(source/'reference_library_pins.json');need('installed reference 9 critical files',len(pins['files'])==9 and raw['reference']['reference_library_sha256']==pins['files'] and (public or before['runtime_critical_file_sha256']==pins['files']))
 replay_max=0.
 def equivalent(x,y):
  nonlocal replay_max
  if isinstance(x,float):
   d=abs(x-y);replay_max=max(replay_max,d);return math.isfinite(d) and d<=1e-12
  if isinstance(x,dict):return x.keys()==y.keys() and all(equivalent(v,y[k]) for k,v in x.items())
  if isinstance(x,list):return len(x)==len(y) and all(equivalent(v,w) for v,w in zip(x,y))
  return x==y
 for k in comp:
  if k!='input_sha256':need('fresh frozen comparator replay '+k,equivalent(comp[k],replay[k]))
 need('replay raw hashes',sorted(comp['input_sha256'].values())==sorted(replay['input_sha256'].values()))
 need('complete science gates',comp['exception'] is None and comp['numerically_eligible'] and comp['completion_count']==2 and not comp['pilot_excluded'] and all(c['passed'] for c in comp['checks']) and {r['id']:r['status'] for r in comp['rows']}=={'P23':'PASS','P24':'PASS'})
 for spec,tc,rc in zip(cfg['cases'],raw['target']['cases'],raw['reference']['cases']):
  cid=spec['id'];K=geom(cfg,spec);M=len(K);derived=[]
  for tag,r in [('target',tc),('reference',rc)]:
   coords=list(zip(r['basis_up'],r['basis_down']));order=sorted(range(len(coords)),key=coords.__getitem__);coords=[coords[i] for i in order];dim=spec['dimension'];nu,nd=spec['Nup'],spec['Ndown'];need(tag+'/'+cid+'/complete basis',len(coords)==len(set(coords))==dim and all(u.bit_count()==nu and d.bit_count()==nd and 0<=u<2**M and 0<=d<2**M for u,d in coords))
   H,C,S=independent_operators(coords,K,spec['U']);V=np.asarray(r['vectors'])[order];E=np.asarray(r['energies']);need(tag+'/'+cid+'/finite full arrays',V.shape==(dim,dim) and E.shape==(dim,) and np.isfinite(V).all() and np.isfinite(E).all() and np.all(np.diff(E)>=0))
   close(tag+'/'+cid+'/independent geometry',K,r['h1'],1e-10)
   for key,m in [('H',H),('C',C),('S2',S)]:close(tag+'/'+cid+'/independent Fock '+key,m,decode(r[key])[np.ix_(order,order)],1e-10)
   close(tag+'/'+cid+'/eigen residual',H@V,V*E,1e-10);close(tag+'/'+cid+'/orthogonality',V.T@V,np.eye(dim),1e-10);close(tag+'/'+cid+'/H S2 commutator',H@S,S@H,1e-10)
   if tag=='target':close(cid+'/physical centering',r['centered_to_physical_shift'],spec['U']*((nu+nd)/2-M/4),1e-10)
   groups=[];st=0
   for j in range(1,dim+1):
    if j==dim or E[j]-E[st]>1e-8:groups.append(np.arange(st,j));st=j
   need(tag+'/'+cid+'/complete8 clusters',len(groups)>=8 and all(E[b[0]]-E[a[-1]]>=1e-6 for a,b in zip(groups,groups[1:])))
   dg=len(groups[0]);need(tag+'/'+cid+'/ground multiplicity',dg==r['ground_dimension']);W0=V[:,groups[0]];X=C@W0;Y=C.T@W0;close(tag+'/'+cid+'/retained raising amplitudes',V.T@X,r['transition_amplitudes'],1e-10);close(tag+'/'+cid+'/retained lowering amplitudes',V.T@Y,r['reverse_transition_amplitudes'],1e-10)
   sub=[]
   for ig,I in enumerate(groups[:8]):
    W=V[:,I];P=W@W.T;spin=np.linalg.eigvalsh(W.T@S@W);allowed=np.array([s*(s+1) for s in range((nu+nd)//2+1)]);close(tag+'/'+cid+f'/spin quantization{ig}',np.min(abs(spin[:,None]-allowed),axis=1),0,1e-8);close(tag+'/'+cid+f'/spin leakage{ig}',S@W,W@(W.T@S@W),1e-10)
    # Projector quadratic form rather than saved per-eigenvector amplitudes.
    sub.append((P,spin,float(np.trace(X.T@P@X)/dg),float(np.trace(Y.T@P@Y)/dg),float(np.mean(E[I])-E[0])))
   total=float(np.sum(X*X)/dg);reverse=float(np.sum(Y*Y)/dg)
   if cid=='DIR_N2_U1':close(tag+'/direction analytic E',E[0],(1-math.sqrt(50))/2,1e-10);close(tag+'/direction raising2',total,2,1e-10);close(tag+'/direction lowering0',reverse,0,1e-10)
   if spec['U']==0:
    levels=np.linalg.eigvalsh(K);close(tag+'/'+cid+'/free ground energy',E[0],sum(levels[:nu])+sum(levels[:nd]),1e-10)
   derived.append((E,sub,total,reverse))
  t,q=derived;close(cid+'/complete external energies',t[0],q[0],1e-8)
  for ig,(x,y) in enumerate(zip(t[1],q[1])):
   for ix,name in enumerate(['projector','spin labels','raising strength','lowering strength','excitation']):close(cid+f'/group{ig}/external '+name,x[ix],y[ix],1e-8)
  metrics.append({'case_id':cid,'dimension':spec['dimension'],'ground_energy':float(t[0][0]),'ground_multiplicity':tc['ground_dimension'],'manifolds':8,'total_raising':t[2],'total_lowering':t[3]})
 need('all56 manifold records',len(comp['manifold_receipts'])==56)
 inputs={str(f.relative_to(p)):sha(f) for f in sorted(p.rglob('*')) if f.is_file() and f!=a.output and f.name!='remote_text_receipts.tar'}
 out={'schema':'PYSCF_TWO_ORBITAL_INDEPENDENT_AUDIT_V1','audit_status':'PASS','integrity':'PASS','catalogue_ids':['P23','P24'],'catalogue_status':{'P23':'PASS','P24':'PASS'},'scientific_completions':2,'pilot_completion_credit':0,'physics_solver_executed':False,'check_count':len(checks),'checks':checks,'complete_case_count':7,'complete_manifold_count':56,'case_metrics':metrics,'source_archive_sha256':ARCHIVE,'source_freeze_sha256':FREEZE,'input_sha256':inputs,'auditor_sha256':sha(__file__),'frozen_comparator_checks':len(comp['checks']),'local_vs_ISAAC_replay_max_float_difference':replay_max,'replay_roundoff_allowance':1e-12,'execution_original_sha256':execution_sha,'execution_evidence_mode':'explicit_public_projection' if a.public_execution_metadata else 'original_raw_public_manifest' if public else 'original_raw_with_remote_transfer_proof'}
 a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({k:out[k] for k in ['audit_status','check_count','scientific_completions','complete_case_count','complete_manifold_count','auditor_sha256']}))
if __name__=='__main__':main()
