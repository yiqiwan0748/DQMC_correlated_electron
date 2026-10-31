#!/usr/bin/env python3
"""Independent QuSpin Hamiltonians/operators, exp_op dynamics, and thermal traces."""
import argparse,contextlib,hashlib,importlib.metadata,io,json,math,os,time,tomllib
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'


def compute(cfg):
    import numpy as np
    from scipy.linalg import eigh,expm
    from quspin.basis import spinful_fermion_basis_1d
    from quspin.operators import hamiltonian,exp_op
    L,U,t=cfg['L'],cfg['U'],cfg['t'];basis=spinful_fermion_basis_1d(L)
    assert L==4 and cfg['boundary']=='pbc'
    def operator(static,hermitian=True,dtype=np.float64):
        with contextlib.redirect_stdout(io.StringIO()):
            return hamiltonian(static,[],basis=basis,dtype=dtype,check_herm=hermitian,check_symm=False,check_pcon=False)
    hops=[]
    for i in range(L):j=(i+1)%L;hops.extend([[-t,i,j],[-t,j,i]])
    Hq=operator([['+-|',hops],['|+-',hops],['n|n',[[U,i,i] for i in range(L)]]]);H=np.asarray(Hq.toarray())
    Nu=[np.asarray(operator([['n|',[[1.,j]]]]).toarray()) for j in range(L)]
    Nd=[np.asarray(operator([['|n',[[1.,j]]]]).toarray()) for j in range(L)]
    D=[np.asarray(operator([['n|n',[[1.,j,j]]]]).toarray()) for j in range(L)]
    NuT=sum(Nu);NdT=sum(Nd);N=NuT+NdT;DT=sum(D)
    nu=np.rint(np.diag(NuT)).astype(int);nd=np.rint(np.diag(NdT)).astype(int)
    diag={'H_hermiticity':float(np.linalg.norm(H-H.T)),'N_commutator':float(np.linalg.norm(H@N-N@H))}
    def real(x):
        x=complex(x);assert abs(x.imag)<1e-9 and math.isfinite(x.real)
        return float(x.real)
    def expect(psi,A):return real(np.vdot(psi,A@psi))
    def block(up,down):
        ids=np.flatnonzero((nu==up)&(nd==down));hh=H[np.ix_(ids,ids)];e,v=eigh(hh)
        return ids,e,v,float(np.linalg.norm(hh@v-v*e))
    # P08: QuSpin exp_op action; independent dense spectral propagation is a gate.
    spec=cfg['P08'];occupied=set(j-1 for j in spec['initial_doublon_sites_one_based'])
    match=np.ones(basis.Ns,dtype=bool)
    for j in range(L):match&=(np.diag(Nu[j])==(j in occupied))&(np.diag(Nd[j])==(j in occupied))
    ids=np.flatnonzero(match);assert len(ids)==1
    initial=np.zeros(basis.Ns,dtype=complex);initial[ids[0]]=1
    times=np.linspace(spec['time_start'],spec['time_stop'],spec['time_count'])
    states=np.asarray(exp_op(Hq,a=-1j,start=times[0],stop=times[-1],num=len(times),endpoint=True,iterate=False).dot(initial))
    assert states.shape==(basis.Ns,len(times))
    e,v=eigh(H);amplitudes=v.T@initial;rows=[];propagation_error=0.
    for index,tm in enumerate(times):
        psi=states[:,index];spectral=v@(np.exp(-1j*tm*e)*amplitudes)
        propagation_error=max(propagation_error,float(np.linalg.norm(psi-spectral)))
        density=[expect(psi,Nu[j]+Nd[j]) for j in range(L)]
        rows.append({'time':float(tm),'norm':real(np.vdot(psi,psi)),'energy':expect(psi,H),
                     'Nup':expect(psi,NuT),'Ndown':expect(psi,NdT),'density':density,
                     'imbalance':sum((1 if j%2 else -1)*density[j] for j in range(L))/(spec['Nup']+spec['Ndown']),
                     'double_occupancy_total':expect(psi,DT)})
    P08={'rows':rows,'exp_action_vs_dense_spectral_error':propagation_error,'initial_energy':expect(initial,H)}
    # P09: complete spin-resolved N+/-1 blocks; poles grouped before comparison.
    spec=cfg['P09'];ids,e,v,residual=block(spec['Nup'],spec['Ndown']);E0=float(e[0])
    psi=np.zeros(basis.Ns);psi[ids]=v[:,0]
    ground={'energy':E0,'norm':real(np.vdot(psi,psi)),'gap_same_sector':float(e[1]-e[0]),
            'residual':residual,'Nup':expect(psi,NuT),'Ndown':expect(psi,NdT)}
    rows=[]
    def clusters(energies,weights):
        order=np.argsort(energies);groups=[];i=0
        while i<len(order):
            indices=[int(order[i])];anchor=energies[order[i]];j=i+1
            while j<len(order) and energies[order[j]]-anchor<=spec['degeneracy_cluster_abs']:indices.append(int(order[j]));j+=1
            groups.append({'omega':float(np.mean(energies[indices])),'weight':float(np.sum(weights[indices])),
                           'multiplicity':len(indices),'span':float(np.ptp(energies[indices]))})
            i=j
        return groups
    mu=spec['physical_mu_frequency'];G=H-mu*N;g0=E0-mu*(spec['Nup']+spec['Ndown'])
    for spin in spec['spin']:
        up0,down0=spec['Nup'],spec['Ndown'];plus=(up0+(spin=='up'),down0+(spin=='down'));minus=(up0-(spin=='up'),down0-(spin=='down'))
        ip,ep,vp,rp=block(*plus);im,em,vm,rm=block(*minus)
        opstr='-|' if spin=='up' else '|-'
        for m in spec['momentum_indices']:
            k=2*np.pi*m/L;c=np.asarray(operator([[opstr,[[np.exp(-1j*k*j)/math.sqrt(L),j] for j in range(L)]]],False,np.complex128).toarray())
            add=c.conj().T@psi;remove=c@psi
            wa=abs(vp.conj().T@add[ip])**2;wr=abs(vm.conj().T@remove[im])**2
            fa=ep-E0-mu;fr=E0-em-mu
            first=float(wa@fa+wr@fr);second=float(wa@(fa*fa)+wr@(fr*fr))
            Xa=(G-g0*np.eye(basis.Ns))@add;Xr=(G-g0*np.eye(basis.Ns))@remove
            direct_first=real(np.vdot(add,Xa)-np.vdot(remove,Xr));direct_second=real(np.vdot(Xa,Xa)+np.vdot(Xr,Xr))
            rows.append({'spin':spin,'momentum_index':m,'addition':clusters(fa,wa),'removal':clusters(fr,wr),
                'addition_dimension':len(ip),'removal_dimension':len(im),'addition_eigen_residual':rp,'removal_eigen_residual':rm,
                'occupation':real(np.vdot(remove,remove)),'addition_norm':real(np.vdot(add,add)),
                'sum_weight':float(sum(wa)+sum(wr)),'first_moment':first,'first_moment_direct':direct_first,
                'first_moment_expected':-2*t*math.cos(k)+U*expect(psi,NdT if spin=='up' else NuT)/L-mu,
                'second_moment':second,'second_moment_direct':direct_second,
                'anticommutator_matrix_error':float(np.linalg.norm(c@c.conj().T+c.conj().T@c-np.eye(basis.Ns)))})
    cross=spec['cross_spin_momentum_pair']
    cu=np.asarray(operator([['-|',[[np.exp(-2j*np.pi*cross[0]*j/L)/math.sqrt(L),j] for j in range(L)]]],False,np.complex128).toarray())
    cd=np.asarray(operator([['|-',[[np.exp(-2j*np.pi*cross[1]*j/L)/math.sqrt(L),j] for j in range(L)]]],False,np.complex128).toarray())
    P09={'ground':ground,'rows':rows,'operator_checks':{
        'cross_spin_annihilation_anticommutator':float(np.linalg.norm(cu@cd+cd@cu)),
        'cross_spin_mixed_anticommutator':float(np.linalg.norm(cu@cd.conj().T+cd.conj().T@cu))}}
    # P10: full-GCE continuum Lehmann sum. Nonzero means are subtracted explicitly.
    spec=cfg['P10'];mu=spec['physical_mu'];G=H-mu*N;e,v=eigh(G);shift=e-e[0];rows=[];thermal=[]
    for beta in spec['beta']:
        p=np.exp(-beta*shift);Z=float(sum(p));p/=Z;rho=(v*p)@v.T
        rho_expm=expm(-beta*(G-e[0]*np.eye(basis.Ns)));rho_expm/=np.trace(rho_expm)
        means=[real(np.trace(rho@(Nu[j]-Nd[j])/2)) for j in range(L)]
        moment=real(np.trace(rho@(N/4-DT/2)))-sum(x*x for x in means)
        thermal.append({'beta':beta,'mu_physical':mu,'mu_centered':mu-U/2,'trace':real(np.trace(rho)),
                        'rho_expm_error':float(np.linalg.norm(rho-rho_expm)),
                        'Nmean':real(np.trace(rho@N)),'double_occupancy_total':real(np.trace(rho@DT)),
                        'local_moment_sum':moment,'spin_mean_abs_max':max(abs(x) for x in means)})
        for m in spec['momentum_indices']:
            q=2*np.pi*m/L;S=sum(np.exp(-1j*q*j)*(Nu[j]-Nd[j])/2 for j in range(L))/math.sqrt(L)
            A=v.conj().T@S@v;weights=abs(A)**2;mean=complex(np.trace(rho@S));offset=abs(mean)**2
            endpoint=real(np.trace(rho@S@S.conj().T))-offset
            for fraction in spec['tau_fractions']:
                tau=beta*fraction;boltz=np.exp(-(beta-tau)*shift[:,None]-tau*shift[None,:])/Z
                value=float(np.sum(boltz*weights))-offset
                rows.append({'beta':beta,'momentum_index':m,'tau_fraction':fraction,'tau':tau,'value':value,
                    'equal_time_direct':endpoint,'spin_mean_abs':abs(mean),
                    'spectral_weight_min':float(weights.min()),'unconnected_value':value+offset})
    P10={'rows':rows,'thermal':thermal}
    return {'P08':P08,'P09':P09,'P10':P10,'diagnostics':diag}


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',default=str(Path(__file__).with_name('manifest.toml')));p.add_argument('--output',required=True)
    args=p.parse_args();path=Path(args.manifest);cfg=tomllib.loads(path.read_text());began=time.monotonic()
    versions={k:importlib.metadata.version(k) for k in ('quspin','numpy','scipy')}
    for k,v in versions.items():assert v==cfg[k+'_version']
    out={'schema':cfg['schema'],'status':'COMPLETED','manifest_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'oracle_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'versions':versions}
    out.update(compute(cfg));out['elapsed_seconds']=time.monotonic()-began
    output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as f:json.dump(out,f,indent=2,allow_nan=False)

if __name__=='__main__':main()
