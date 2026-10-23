#!/usr/bin/env python3
"""Independent proof/reproduction of P16; does not change a solver or old result."""
import argparse,contextlib,hashlib,io,json,itertools,os,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from quspin.basis import spinful_fermion_basis_1d
from quspin.operators import hamiltonian

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def op(b,terms):
    with contextlib.redirect_stdout(io.StringIO()):
        h=hamiltonian(terms,[],basis=b,dtype=np.float64,check_herm=False,check_symm=False,check_pcon=False)
    return np.asarray(h.toarray())
def phase(up,dn,A,B):
    within=sum(((up>>i)&1)*((up>>j)&1)+((dn>>i)&1)*((dn>>j)&1) for i in B for j in A if i<j)
    between=sum((up>>i)&1 for i in B)*sum((dn>>i)&1 for i in A)
    return -1. if (within+between)%2 else 1.
def explicit_phase(up,dn,A,B,L):
    modes=A+[L+i for i in A]+B+[L+i for i in B]
    present=[i for i in modes if ((up|(dn<<L))>>i)&1]
    return -1. if sum(i>j for k,i in enumerate(present) for j in present[k+1:])%2 else 1.
def compact(mask,A):return sum(((mask>>i)&1)<<k for k,i in enumerate(A))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.time();L=4;U=4.;checks=0
    for bits in range(1<<L):
        A=[i for i in range(L) if bits>>i&1];B=[i for i in range(L) if i not in A]
        for up in range(1<<L):
            for dn in range(1<<L):
                assert phase(up,dn,A,B)==explicit_phase(up,dn,A,B,L)
                checks+=1
    b=spinful_fermion_basis_1d(L);hops=[[-1.,i,(i+1)%L] for i in range(L)]+[[-1.,(i+1)%L,i] for i in range(L)]
    H=op(b,[['+-|',hops],['|+-',hops],['n|n',[[U,i,i] for i in range(L)]]])
    nu=np.array([np.diag(op(b,[['n|',[[1.,i]]]])) for i in range(L)]);nd=np.array([np.diag(op(b,[['|n',[[1.,i]]]])) for i in range(L)])
    idx=np.where((nu.sum(0)==2)&(nd.sum(0)==2))[0];e,v=eigh(H[np.ix_(idx,idx)]);psi=np.zeros(b.Ns);psi[idx]=v[:,0]
    rows=[]
    for bits in range(1<<L):
        A=[i for i in range(L) if bits>>i&1];B=[i for i in range(L) if i not in A];na=len(A);nb=len(B)
        X=np.zeros((1<<(2*na),1<<(2*nb)));old=np.zeros_like(X)
        for k in idx:
            up=sum(int(round(nu[i,k]))<<i for i in range(L));dn=sum(int(round(nd[i,k]))<<i for i in range(L))
            ua,ub,da,db=compact(up,A),compact(up,B),compact(dn,A),compact(dn,B)
            X[ua|(da<<na),ub|(db<<nb)]+=phase(up,dn,A,B)*psi[k]
            old[ua|(da<<na),ub|(db<<nb)]+=(-1. if (ub.bit_count()*da.bit_count())%2 else 1.)*psi[k]
        pc=float(np.sum((X@X.T)**2));po=float(np.sum((old@old.T)**2))
        if not A or not B:pq=1.
        else:
            with contextlib.redirect_stdout(io.StringIO()):rho=b.partial_trace(psi,sub_sys_A=(tuple(A),tuple(A)),return_rdm='A',enforce_pure=True)
            rho=np.squeeze(np.asarray(rho));pq=float(np.trace(rho@rho))
        rows.append({'sites':A,'old_purity':po,'corrected_purity':pc,'quspin_purity':pq,'corrected_absolute_error':abs(pc-pq)})
    out={'purpose':'P16 diagnostic and candidate correction, not an additional catalogue idea','job_id':os.getenv('SLURM_JOB_ID','LOCAL'),'source_sha256':sha(__file__),'formula_permutation_checks':checks,'ground_energy':float(e[0]),'ground_gap':float(e[1]-e[0]),'ground_residual':float(np.linalg.norm(H@psi-e[0]*psi)),'all_subset_purities':rows,'maximum_corrected_quspin_error':max(x['corrected_absolute_error'] for x in rows),'elapsed_seconds':time.time()-start,'status':'PASS' if max(x['corrected_absolute_error'] for x in rows)<=1e-10 else 'FAIL','quspin_source':'https://quspin.github.io/QuSpin/modules/quspin/basis/basis_1d/fermion.html'}
    with a.output.open('x') as f:json.dump(out,f,indent=2)
    print(json.dumps({k:v for k,v in out.items() if k!='all_subset_purities'}));sys.exit(out['status']!='PASS')
if __name__=='__main__':main()
