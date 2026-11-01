#!/usr/bin/env python3
"""Independent postprocessing of retained X17 states; no native/QuSpin solver call."""
import argparse,hashlib,itertools,json,math,tomllib,traceback
from pathlib import Path
import numpy as np
FREEZE='171c58736e7f0e59443b4429848827e35e81e0fa44c84ac2eb53dddd15e8fbcb'
ARCHIVE='5dfe6b1718a9af77fbfde134d09194c7594c1a144ce0d3bb6314fabaa3b13958'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return tomllib.loads(p.read_text()) if p.suffix=='.toml' else json.loads(p.read_text())
def hop(f,p,q):
 if not f>>q&1:return None
 sign=(-1)**((f&((1<<q)-1)).bit_count());f^=1<<q
 if f>>p&1:return None
 return f|(1<<p),sign*(-1)**((f&((1<<p)-1)).bit_count())
def geometry(cfg,t4):
 K=np.zeros((8,8));sites=[(x,y,o) for y in range(2) for x in range(2) for o in range(2)]
 for i,(x,y,o) in enumerate(sites):
  for j,(a,b,c) in enumerate(sites):
   dx,dy=a-x,b-y
   if o==c and abs(dx)==1 and dy==0:K[i,j]=-cfg['t1' if o==0 else 't2']
   if o==c and dx==0 and abs(dy)==1:K[i,j]=-cfg['t2' if o==0 else 't1']
   if abs(dx)==abs(dy)==1:K[i,j]=-cfg['t3'] if o==c else t4*dx*dy
 return K

def hamiltonian(coords,K,U):
 ix={u|(d<<8):i for i,(u,d) in enumerate(coords)};H=np.zeros((len(coords),len(coords)))
 for j,(u,d) in enumerate(coords):
  H[j,j]=U*((u&d).bit_count()-(u.bit_count()+d.bit_count())/2+2)
  for sigma in (0,1):
   for i,k in zip(*np.nonzero(K)):
    r=hop(u|(d<<8),int(i)+8*sigma,int(k)+8*sigma)
    if r:H[ix[r[0]],j]+=K[i,k]*r[1]
 return H

def decode(row):
 a=np.zeros(row['shape']);assert len(set(zip(row['row'],row['col'])))==len(row['value']);a[row['row'],row['col']]=row['value'];return a

def explicit_partition(psi,coords,A,B):
 order=A+[8+i for i in A]+B+[8+i for i in B];X=np.zeros((256,256))
 for coeff,(u,d) in zip(psi,coords):
  f=u|(d<<8);occupied=[i for i in order if f>>i&1];parity=sum(a>b for j,a in enumerate(occupied) for b in occupied[j+1:]);sign=(-1)**parity
  a=sum(((u>>site)&1)<<j for j,site in enumerate(A))+sum(((d>>site)&1)<<(j+4) for j,site in enumerate(A));b=sum(((u>>site)&1)<<j for j,site in enumerate(B))+sum(((d>>site)&1)<<(j+4) for j,site in enumerate(B));X[a,b]+=sign*coeff
 return X@X.T,X.T@X

def observables(psi,coords,A,B):
 focks=[u|(d<<8) for u,d in coords];where={f:i for i,f in enumerate(focks)};up=np.array([[(u>>i)&1 for u,d in coords]for i in range(8)]);dn=np.array([[(d>>i)&1 for u,d in coords]for i in range(8)]);n=up+dn;z=(up-dn)/2;p=psi*psi;nm=n@p;zm=z@p;nn=np.zeros((4,4));sd=np.zeros((4,4));zz=np.zeros((4,4));Sp={}
 for k,f in enumerate(focks):
  for i in range(8):
   r=hop(f,i,i+8)
   if r:Sp[r[0]]=Sp.get(r[0],0.)+r[1]*psi[k]
 for i,a in enumerate(A):
  for j,b in enumerate(B):
   nn[i,j]=(n[a]*n[b])@p-nm[a]*nm[b];zz[i,j]=(z[a]*z[b])@p-zm[a]*zm[b];sd[i,j]=zz[i,j]
   for k,f in enumerate(focks):
    for bp,bq,ap,aq in [(b+8,b,a,a+8),(b,b+8,a+8,a)]:
     r=hop(f,bp,bq)
     if r:
      s=hop(r[0],ap,aq)
      if s:sd[i,j]+=.5*psi[where[s[0]]]*psi[k]*r[1]*s[1]
 return dict(density_up=up@p,density_down=dn@p,density_total=nm,Sz_mean=zm,cross_nn_connected=nn,cross_SzSz_connected=zz,cross_Sdot_connected=sd,spin_squared=sum(x*x for x in Sp.values()))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--receipt-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();p=a.receipt_root;assert not a.output.exists();checks=[];metrics=[];out={}
 def need(n,b):checks.append(dict(name=n,passed=bool(b)));assert b,n
 def close(n,x,y,t):
  x=np.asarray(x);y=np.asarray(y);assert x.shape==y.shape,n;err=float(np.max(abs(x-y),initial=0));checks.append(dict(name=n,error=err,limit=t,passed=math.isfinite(err) and err<=t));assert checks[-1]['passed'],n;return err
 try:
  source=p/'source';cfg=load(source/'manifest.toml');g=cfg['gates'];freeze=load(source/'source_freeze.json');need('exact original source freeze',sha(source/'source_freeze.json')==FREEZE and freeze['archive_sha256']==ARCHIVE)
  for n,h in freeze['files'].items():need('source '+n,sha(source/n)==h)
  before=load(p/'transfer_before.json');after=load(p/'transfer_after.json');aft={r['local']:r for r in after['files']};need('remote-after stable',after['status']=='PASS')
  for r in before['files']:
   need('input hash '+r['local'],sha(p/r['local'])==r['sha256'] and (p/r['local']).stat().st_size==r['bytes']);v=aft[r['local']];need('remote unchanged '+r['local'],v['sha256']==r['sha256'] and v['mtime_ns']==r['before']['mtime_ns'] and v['inode']==r['before']['inode'])
  target=load(p/'results/target.toml');oracle=load(p/'results/reference.json');saved=load(p/'results/comparison.json');fresh=load(p/'frozen_pilot_audit/comparison_replay.json');resources=load(p/'frozen_pilot_audit/admission.json');execution=load(p/'results/execution.json')
  need('successful pilot processes',execution['mode']=='pilot' and execution['status']=='PROCESSES_COMPLETE' and [x['name']for x in execution['processes']]==['target','reference','comparator'] and all(x['exit_code']==0 and not x['timed_out']for x in execution['processes']))
  need('pilot0credit',saved['scientific_completions']==fresh['scientific_completions']==0 and saved['catalogue_status']==fresh['catalogue_status']=='NOT_RUN' and saved['pilot_excluded'] and fresh['pilot_excluded'])
  need('frozen arithmetic replay all gates',saved['status']==fresh['status']=='PASS' and all(x['passed']for c in fresh['cases']for x in c['gates']))
  library=load(source/'reference_library_pins.json');need('all72 reference pins',len(library['files'])==72 and oracle['reference_library_sha256']==before['runtime_critical_file_sha256']==library['files'])
  specs=[s for s in cfg['cases'] if s['id'] in cfg['pilot_case_ids']];need('all pilot fixtures',[s['id']for s in target['cases']]==[s['id']for s in oracle['cases']]==cfg['pilot_case_ids'])
  for spec,t,r in zip(specs,target['cases'],oracle['cases']):
   cid=spec['id'];K=geometry(cfg,spec['t4']);close(cid+'/independent geometry target',K,t['h1'],g['hamiltonian_abs']);close(cid+'/independent geometry oracle',K,r['h1'],g['hamiltonian_abs']);need(cid+'/both producers ready',t['status']=='TARGET_READY' and r['status']=='ORACLE_READY');Eall=[];Eblocks=[];budget=[];primary=None
   need(cid+'/all5 sectors',[(s['Nup'],s['Ndown'])for s in t['sectors']]==[(s['Nup'],s['Ndown'])for s in r['sectors']]==[tuple(x)for x in cfg['spin_sectors']])
   for st,sr,d in zip(t['sectors'],r['sectors'],cfg['sector_dimensions']):
    nu,nd=st['Nup'],st['Ndown'];coords=[(u,v)for u in st['basis_up_masks']for v in st['basis_down_masks']];need(cid+f'/basis{nu}',len(coords)==len(set(coords))==d and all(u.bit_count()==nu and v.bit_count()==nd and 0<=u<256 and 0<=v<256 for u,v in coords));H=hamiltonian(coords,K,spec['U']);close(cid+f'/independent H{nu}',H,decode(st['H']),g['hamiltonian_abs']);refcoords=list(zip(sr['basis_up_per_row'],sr['basis_down_per_row']));where={x:i for i,x in enumerate(refcoords)};perm=[where[x]for x in coords];close(cid+f'/independent oracle H{nu}',H,decode(sr['H'])[np.ix_(perm,perm)],g['hamiltonian_abs'])
    E,V=np.linalg.eigh(H);close(cid+f'/complete native spectrum{nu}',E,st['energies'],g['energy_abs']);close(cid+f'/complete oracle spectrum{nu}',E,sr['energies'],g['energy_abs']);Eall.extend(E);Eblocks.append(E);scale=max(1.,np.linalg.norm(H,np.inf));res=max(np.linalg.norm(H@V-V*E,axis=0));orth=np.max(abs(V.T@V-np.eye(d)));budget.append(2*np.sqrt(d)*res+2*d*orth*scale+g['gap_roundoff_multiplier']*np.finfo(float).eps*d*scale)
    psi=np.asarray(st['native_vector']);close(cid+f'/native norm{nu}',np.dot(psi,psi),1.,g['rdm_trace_abs']);close(cid+f'/native physical residual{nu}',np.linalg.norm(H@psi-st['native_energy']*psi)/scale,0.,g['residual_scaled']);ng=np.count_nonzero(E-E[0]<=g['ground_cluster_diameter']);fidelity=float(np.sum(abs(V[:,:ng].T@psi)**2));need(cid+f'/native ground fidelity{nu}',fidelity>=g['native_ground_fidelity_min'] and st['native_converged'] and st['native_method']=='lanczos')
    if nu==2:primary=(coords,H,E,V,psi,fidelity)
   coords,H,E,V,psi,fid=primary;ordered=np.sort(Eall);need(cid+'/unique global ground',E[1]-E[0]>g['global_gap_min'] and abs(E[0]-ordered[0])<=g['energy_abs'] and ordered[1]-ordered[0]>g['global_gap_min']);lower=min(float(e)-budget[j]for j,block in enumerate(Eblocks)for k,e in enumerate(block)if not(j==2 and k==0))-(E[0]+budget[2]);need(cid+'/independent conditional gap',lower>g['global_gap_min'])
   A=cfg['orbital_A_zero_based'];B=cfg['orbital_B_zero_based'];ra,rb=explicit_partition(psi,coords,A,B);obs=observables(psi,coords,A,B)
   for part,rho in [('A',ra),('B',rb)]:
    vals=np.linalg.eigvalsh(rho);close(cid+'/'+part+'/explicit fermion RDM',rho,t['partition']['rdm_'+part],g['rdm_hermiticity_abs']);close(cid+'/'+part+'/native full256 levels',vals,t['partition']['levels_'+part],g['purity_abs']);l1=float(np.sum(abs(vals-np.asarray(r['partition']['levels_'+part]))));need(cid+'/'+part+'/external full256 L1',l1<=g['rdm_spectrum_l1']);pur=float(np.sum(rho*rho));s2=float(-np.log(pur));close(cid+'/'+part+'/actual native entropy',s2,t['partition']['native_S2_'+part],g['entropy_abs']);close(cid+'/'+part+'/external entropy',s2,r['partition']['S2_'+part],g['entropy_abs'])
   for key in ['density_up','density_down','density_total','Sz_mean','cross_nn_connected','cross_SzSz_connected','cross_Sdot_connected','spin_squared']:
    close(cid+'/independent native '+key,obs[key],t['observables'][key],g['connected_observable_abs']);close(cid+'/independent external '+key,obs[key],r['observables'][key],g['connected_observable_abs'])
   need(cid+'/independent singlet',obs['spin_squared']<=g['spin_squared_abs']);metrics.append(dict(case_id=cid,complete_spin_sectors=5,balanced_native_ground_fidelity=fid,independent_global_gap=float(ordered[1]-ordered[0]),independent_conditional_gap_lower=float(lower),explicit_permutation_amplitudes=len(psi),orbital_levels_per_cut=256,entropy_A=float(-np.log(np.sum(ra*ra))),max_abs_connected_density=float(np.max(abs(obs['cross_nn_connected']))),max_abs_connected_Sdot=float(np.max(abs(obs['cross_Sdot_connected'])))))
  def rss(f):
   import re
   return int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',f.read_text())[1])
  peaks=[rss(p/'results'/f'{q}.resources.txt')for q in ['target','reference','comparator']];secs=sum(q['elapsed_seconds']for q in execution['processes']if q['name']!='comparator');proj=1.5*6*secs+60;memory=max(peaks)+15*((p/'results/target.toml').stat().st_size+(p/'results/reference.json').stat().st_size)/1024;close('unchanged projection seconds',proj,resources['projected_full_seconds'],1e-8);close('unchanged projection RSS',memory,resources['projected_full_peak_rss_KiB'],1e-8);need('engineering eligibility',resources['status']=='PASS' and max(peaks)<=2800000 and memory<=3000000 and proj<=1600)
  out.update(audit_status='PASS',integrity='PASS',catalogue_status='NOT_RUN',scientific_completions=0,source_archive_sha256=ARCHIVE,source_freeze_sha256=FREEZE,case_metrics=metrics,projected_full_seconds=proj,projected_full_peak_rss_KiB=memory,resource_admission_only=True,root_production_admitted=False)
 except Exception as ex:out.update(audit_status='UNRESOLVED',catalogue_status='NOT_RUN',scientific_completions=0,reason=type(ex).__name__+': '+str(ex),traceback=traceback.format_exc(),case_metrics=metrics)
 out.update(schema='X17_PILOT_INDEPENDENT_RAW_AUDIT_V1',checks=checks,check_count=len(checks),scientific_producer_rerun=False,postprocessing_eigendecomposition=True,input_sha256={str(f.relative_to(p)):sha(f)for f in sorted(p.rglob('*'))if f.is_file() and f!=a.output and f.name!='remote_text_receipts.tar'},auditor_sha256=sha(__file__))
 a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({k:out[k]for k in ['audit_status','catalogue_status','scientific_completions','check_count']}))
if __name__=='__main__':main()
