#!/usr/bin/env python3
"""Independent full-trace QuSpin oracle for three ensemble/response questions."""
import argparse
import contextlib
import hashlib
import importlib.metadata
import io
import json
import math
import os
from pathlib import Path
import time
import tomllib

for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'


def run(cfg):
    import numpy as np
    from scipy.linalg import eigh
    from scipy.special import logsumexp
    from quspin.basis import spinful_fermion_basis_1d
    from quspin.operators import hamiltonian
    L,U,t,nA=cfg['L'],cfg['U'],cfg['t'],cfg['nA']
    assert L==4 and cfg['boundary']=='pbc'
    basis=spinful_fermion_basis_1d(L)
    def matrix(static):
        with contextlib.redirect_stdout(io.StringIO()):
            h=hamiltonian(static,[],basis=basis,dtype=np.float64,
                          check_herm=True,check_symm=False,check_pcon=False)
        return np.asarray(h.toarray(),dtype=float)
    hops=[]
    for i in range(L):
        j=(i+1)%L;hops.extend([[-t,i,j],[-t,j,i]])
    K=matrix([['+-|',hops],['|+-',hops]])
    D=matrix([['n|n',[[1.,i,i] for i in range(L)]]])
    Nu=matrix([['n|',[[1.,i] for i in range(L)]]])
    Nd=matrix([['|n',[[1.,i] for i in range(L)]]])
    N=Nu+Nd;H=K+U*D
    n=np.rint(np.diag(N)).astype(int)
    up=np.rint(np.diag(Nu)).astype(int);down=np.rint(np.diag(Nd)).astype(int)
    assert np.linalg.norm(H@N-N@H)<1e-12
    blocks={}
    for number in range(2*L+1):
        indices=np.flatnonzero(n==number);e,v=eigh(H[np.ix_(indices,indices)])
        blocks[number]=(indices,e,v)
    levels=np.concatenate([block[1] for block in blocks.values()])
    charges=np.concatenate([np.full(len(block[1]),number) for number,block in blocks.items()])
    diagnostics={'maximum_trace_error':0.,'maximum_rdm_trace_error':0.,
                 'hermiticity_error':float(np.linalg.norm(H-H.T)),
                 'N_commutator_error':float(np.linalg.norm(H@N-N@H))}
    def real(x):
        x=complex(x);assert abs(x.imag)<1e-10 and math.isfinite(x.real)
        return float(x.real)
    def expect(rho,A):return real(np.trace(rho@A))
    def density(beta,mu):
        e,v=eigh(H-mu*N)
        p=np.exp(-beta*(e-e[0]));p/=p.sum()
        rho=(v*p)@v.T
        diagnostics['maximum_trace_error']=max(diagnostics['maximum_trace_error'],abs(real(np.trace(rho))-1))
        return rho,e,v,p
    def purity(rho):
        r=np.asarray(basis.partial_trace(rho,sub_sys_A=(list(range(nA)),list(range(nA))),return_rdm='A',enforce_pure=False))
        if r.shape==(1,4**nA,4**nA):r=r[0]
        assert r.shape==(4**nA,4**nA)
        diagnostics['maximum_rdm_trace_error']=max(diagnostics['maximum_rdm_trace_error'],abs(real(np.trace(r))-1))
        p=real(np.trace(r@r));return p,-math.log(p)
    def number_weight(rho,Ntarget):return real(np.diag(rho)[n==Ntarget].sum())
    def spin_weight(rho,Nu0,Nd0):return real(np.diag(rho)[(up==Nu0)&(down==Nd0)].sum())
    def moments(rho,mu):
        G=H-mu*N;E=expect(rho,H);Nm=expect(rho,N);Gm=expect(rho,G)
        E2=expect(rho,H@H);N2=expect(rho,N@N);G2=expect(rho,G@G)
        EN=expect(rho,H@N);EG=expect(rho,H@G)
        return {'Ephysical':E,'Nmean':Nm,'Gphysical':Gm,'E2':E2,'N2':N2,'G2':G2,'EN':EN,'EG':EG,
                'var_E':E2-E*E,'var_N':N2-Nm*Nm,'var_G':G2-Gm*Gm,
                'cov_EN':EN-E*Nm,'cov_EG':EG-E*Gm}
    # X02: unique GCE ground state, independently checked against fixed-total-N minima.
    spec=cfg['X02'];mu=spec['mu_physical'];N0=spec['N_ground']
    _,e,v,_=density(1.,mu);psi=v[:,0];rho0=np.outer(psi,psi)
    p0,s0=purity(rho0);E0=expect(rho0,H)
    ground={'Ephysical':E0,'S2':s0,'purity':p0,'Nmean':expect(rho0,N),
            'Nup_mean':expect(rho0,Nu),'Ndown_mean':expect(rho0,Nd),
            'global_grand_gap':float(e[1]-e[0]),
            'ground_residual':float(np.linalg.norm((H-mu*N)@psi-e[0]*psi)),
            'mu_minus':float(blocks[N0][1][0]-blocks[N0-1][1][0]),
            'mu_plus':float(blocks[N0+1][1][0]-blocks[N0][1][0])}
    rows=[]
    for beta in spec['beta']:
        rho,_,_,p=density(beta,mu);pur,s2=purity(rho);epsilon=float(p[1:].sum())
        B=spec['purity_lipschitz_factor']*epsilon;valid=B<p0
        bound=-math.log1p(-B/p0) if valid else -1.
        wN=number_weight(rho,N0);delta=s2-s0
        qualifies=bool(valid and bound<=spec['entropy_bias_target'] and max(0.,1-wN)<=spec['particle_sector_leakage_target'] and abs(delta)<=spec['entropy_bias_target'])
        rows.append({'beta':beta,'mu_physical':mu,'mu_centered':mu-U/2,'S2':s2,'purity':pur,
                     'delta_S2':delta,'Nmean':expect(rho,N),'weight_N_ground':wN,
                     'weight_fixed_Nup_Ndown':spin_weight(rho,spec['Nup_ground'],spec['Ndown_ground']),
                     'ground_weight':float(p[0]),'excited_weight':epsilon,
                     'entropy_bound_valid':valid,'entropy_bias_bound':bound,'qualifies':qualifies})
    X02={'ground':ground,'rows':rows}
    # X03: total-N canonical density includes every allowed spin sector.
    spec=cfg['X03'];Nc=spec['N_total'];indices,ec,vc=blocks[Nc];rows=[]
    for beta in spec['beta']:
        pc=np.exp(-beta*(ec-ec[0]));pc/=pc.sum()
        canonical=np.zeros_like(H);canonical[np.ix_(indices,indices)]=(vc*pc)@vc.T
        lo,hi=spec['mu_physical_bracket']
        def mean(mu):
            z=-beta*(levels-mu*charges);p=np.exp(z-logsumexp(z))
            return float(p@charges),float(beta*((p@(charges*charges))-(p@charges)**2))
        assert mean(lo)[0]<Nc<mean(hi)[0]
        for iteration in range(spec['max_bisections']):
            mid=(lo+hi)/2
            if mean(mid)[0]<Nc:lo=mid
            else:hi=mid
            if hi-lo<=spec['mu_bracket_width']:break
        mu=(lo+hi)/2;gce,_,_,_=density(beta,mu)
        pcA,scA=purity(canonical);pgA,sgA=purity(gce)
        cE=expect(canonical,H)/L;gE=expect(gce,H)/L
        rows.append({'beta':beta,'mu_physical':mu,'mu_centered':mu-U/2,'Nmean_GCE':expect(gce,N),
                     'Nmean_canonical':expect(canonical,N),'Nvariance_canonical':expect(canonical,N@N)-Nc*Nc,
                     'compressibility':mean(mu)[1],'mean_N_error':abs(expect(gce,N)-Nc),'bisections':iteration+1,
                     'mu_bracket_width':hi-lo,'canonical_dimension':len(indices),
                     'fixed_Nup2_Ndown1_dimension':int(np.sum((up==2)&(down==1))),
                     'canonical_weight_Nup2_Ndown1':spin_weight(canonical,2,1),
                     'Ecanonical_per_site':cE,'EGCE_per_site':gE,'delta_energy_per_site':gE-cE,
                     'purity_canonical':pcA,'purity_GCE':pgA,'delta_purity':pgA-pcA,
                     'S2_canonical':scA,'S2_GCE':sgA})
    X03={'rows':rows,'ensemble':'fixed TOTAL N; all Sz sectors included'}
    # P07: differentiate at fixed PHYSICAL mu, never retune it with temperature.
    spec=cfg['P07'];mu=spec['mu_physical'];rows=[]
    for beta in spec['beta']:
        rho,_,_,_=density(beta,mu);m=moments(rho,mu);T=1/beta;h=T*spec['finite_difference_relative_T_step']
        def energy_at(temp):return expect(density(1/temp,mu)[0],H)
        fd=(energy_at(T-2*h)-8*energy_at(T-h)+8*energy_at(T+h)-energy_at(T+2*h))/(12*h)
        response=beta*beta*m['cov_EG']
        rows.append(dict(m,beta=beta,mu_physical=mu,mu_centered=mu-U/2,
                    C_mu=response,beta2_var_grand=beta*beta*m['var_G'],C_mu_finite_difference=fd,
                    finite_difference_error=abs(response-fd),
                    covariance_identity_error=abs(m['cov_EG']-(m['var_E']-mu*m['cov_EN']))))
    return {'X02':X02,'X03':X03,'P07':{'rows':rows},'diagnostics':diagnostics}


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',default=str(Path(__file__).with_name('manifest.toml')));p.add_argument('--output',required=True)
    args=p.parse_args();cfg=tomllib.loads(Path(args.manifest).read_text());start=time.monotonic()
    versions={x:importlib.metadata.version(x) for x in ('quspin','numpy','scipy')}
    for package in ('quspin','numpy','scipy'):
        assert versions[package]==cfg[package+'_version'],f'{package}: unpinned version'
    result={'schema':cfg['schema'],'status':'COMPLETED','manifest_sha256':hashlib.sha256(Path(args.manifest).read_bytes()).hexdigest(),
            'oracle_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'backend_versions':versions}
    result.update(run(cfg));result['elapsed_seconds']=time.monotonic()-start
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)

if __name__=='__main__':main()
