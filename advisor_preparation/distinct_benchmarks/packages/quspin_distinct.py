#!/usr/bin/env python3
"""Independent QuSpin observables; imports no user's solver or target adapter."""
import argparse,contextlib,hashlib,importlib.metadata,io,json,os,platform,time,tomllib
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,eigvalsh
from quspin.basis import spinful_fermion_basis_1d
from quspin.operators import hamiltonian

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def op(b,terms):
 with contextlib.redirect_stdout(io.StringIO()):h=hamiltonian(terms,[],basis=b,dtype=np.float64,check_herm=False,check_symm=False,check_pcon=False)
 return np.asarray(h.toarray())
def flatten(a):return np.asarray(a).reshape(-1).real.tolist()
def spectrum(H,S2,window):
 e,v=eigh(H);ans=[];j=0
 while j<len(e):
  k=j+1
  while k<len(e) and abs(e[k]-e[j])<1e-8:k+=1
  if e[j]-e[0]<=window+1e-9:
   s=eigvalsh(v[:,j:k].T@S2@v[:,j:k]);ans.extend(z for x in s for z in [float(np.mean(e[j:k])-e[0]),float(x)])
  j=k
 return ans

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--manifest',type=Path,default=Path(__file__).with_name('manifest.toml'));ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();a.output_dir.mkdir(exist_ok=False,parents=True);start=time.time();m=tomllib.loads(a.manifest.read_text());L=m['L'];U=m['U'];b=spinful_fermion_basis_1d(L)
 hop=[[float(-m['t']),i,j] for i in range(L) for j in [(i+1)%L]]+[[float(-m['t']),j,i] for i in range(L) for j in [(i+1)%L]]
 H=op(b,[['+-|',hop],['|+-',hop],['n|n',[[float(U),i,i] for i in range(L)]]])
 nu=np.array([np.diag(op(b,[['n|',[[1.,i]]]])) for i in range(L)]);nd=np.array([np.diag(op(b,[['|n',[[1.,i]]]])) for i in range(L)]);N=(nu+nd).sum(axis=0);Sz=(nu-nd).sum(axis=0)/2
 idx=np.where((abs(nu.sum(axis=0)-m['Nup'])<.1)&(abs(nd.sum(axis=0)-m['Ndown'])<.1))[0];e,v=eigh(H[np.ix_(idx,idx)]);assert e[1]-e[0]>m['ground_degeneracy_gap_min'];psi=np.zeros(b.Ns);psi[idx]=v[:,0];assert np.linalg.norm(H@psi-e[0]*psi)<1e-9
 def expect(O):return float(np.real(psi@O@psi))
 def rdm(sites):
  with contextlib.redirect_stdout(io.StringIO()):rho=b.partial_trace(psi,sub_sys_A=(tuple(sites),tuple(sites)),return_rdm='A',enforce_pure=True)
  rho=np.squeeze(np.asarray(rho));assert rho.ndim==2 and abs(np.trace(rho)-1)<1e-10;return rho
 def projections(C):
  w,V=eigh(C);out=w.tolist();j=0
  while j<L:
   k=j+1
   while k<L and abs(w[k]-w[j])<1e-7:k+=1
   out.extend(flatten(V[:,j:k]@V[:,j:k].T));j=k
  return out
 Sp=op(b,[['+|-',[[1.,i,i] for i in range(L)]]]);Sm=Sp.T;S2=np.diag(Sz**2)+(Sp@Sm+Sm@Sp)/2
 assert np.linalg.norm(Sp@Sm-Sm@Sp-2*np.diag(Sz))<1e-10
 assert np.linalg.norm(H@S2-S2@H)<1e-9
 prob=psi**2;n=nu+nd;z=(nu-nd)/2;Czz=(z*prob)@z.T;density=n@prob;Cnn=(n*prob)@n.T-np.outer(density,density)
 def structure(C):return [float(np.real(sum(np.exp(2j*np.pi*q*(i-j)/L)*C[i,j] for i in range(L) for j in range(L))/L)) for q in range(L)]
 Cspin=[np.array([[expect(op(b,[[string,[[1.,i,j]]]])) for j in range(L)] for i in range(L)]) for string in ['+-|','|+-']]
 def pair_matrix(kind):
  operators=[]
  for i in range(L):
   j=(i+1)%L;couplings=[[-1.,i,i]] if kind=='onsite' else [[-1/np.sqrt(2),i,j],[-1/np.sqrt(2),j,i]]
   A=op(b,[['-|-',couplings]]);assert np.linalg.norm(Sp@A-A@Sp)<1e-10;operators.append(A@psi)
  return np.array([[float(x@y) for y in operators] for x in operators])
 values={};errors={}
 for c in m['cases']:
  cid=c['id']
  if not c['ready']:continue
  try:
   if cid=='P01':
    p=np.maximum(eigvalsh(rdm([0,1])),0.);traces=[float(sum(p**q)) for q in [2,3,4]];vals=traces+[float(np.log(x)/(1-q)) for q,x in zip([2,3,4],traces)]
   elif cid=='P02':
    ix=np.where(abs(N-L)<.1)[0];vals=spectrum(H[np.ix_(ix,ix)],S2[np.ix_(ix,ix)],m['spin_window_above_ground'])
   elif cid=='P03':vals=sum([flatten(C)+projections(C) for C in Cspin],[])
   elif cid=='P04':sq=structure(Czz);vals=flatten(Czz)+sq+[float(np.argmax(sq))]
   elif cid=='P05':vals=flatten(Cnn)+structure(Cnn)
   elif cid=='P06':vals=flatten(pair_matrix('onsite'))+flatten(pair_matrix('bond'))
   elif cid=='P16':
    pur=[float(np.trace(rdm(s)@rdm(s))) for s in [[0],[2],[0,2]]];vals=pur+[-np.log(pur[0])-np.log(pur[1])+np.log(pur[2])]
   values[cid]=vals
  except Exception as ex:errors[cid]=f'{type(ex).__name__}: {ex}'
 out=[]
 for cid,vals in values.items():out += [f'[{cid}]','status="ORACLE_READY"','values='+json.dumps(vals),'']
 for cid,error in errors.items():out += [f'[{cid}]','status="UNRESOLVED"','reason='+json.dumps(error),'']
 (a.output_dir/'oracle.toml').write_text('\n'.join(out));receipt={'status':'EXPORTED','manifest_sha256':sha(a.manifest),'source_sha256':sha(__file__),'package_versions':{n:importlib.metadata.version(n) for n in ['quspin','numpy','scipy','quspin-extensions']},'python':platform.python_version(),'job_id':os.getenv('SLURM_JOB_ID','LOCAL'),'elapsed_seconds':time.time()-start,'ground_energy_uncentered':float(e[0]),'ground_residual':float(np.linalg.norm(H@psi-e[0]*psi)),'ground_gap':float(e[1]-e[0]),'sources':['https://quspin.github.io/QuSpin/generated/quspin.basis.spinful_fermion_basis_1d.html','https://quspin.github.io/QuSpin/examples/example6.html'],'errors':errors}
 (a.output_dir/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'rows':list(values),'errors':errors,'seconds':receipt['elapsed_seconds']}))
if __name__=='__main__':main()
