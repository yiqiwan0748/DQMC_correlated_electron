#!/usr/bin/env python3
"""Independent QuSpin ground projectors and prescribed degenerate ensembles."""
import argparse,contextlib,hashlib,importlib.metadata,io,json,math,os,time,tomllib
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'

def compute(cfg):
    import numpy as np
    from scipy.linalg import eigh
    from quspin.basis import spinful_fermion_basis_1d
    from quspin.operators import hamiltonian
    def model(s):
        L=s['L'];basis=spinful_fermion_basis_1d(L)
        def op(static,hermitian=True,dtype=np.float64):
            with contextlib.redirect_stdout(io.StringIO()):
                return np.asarray(hamiltonian(static,[],basis=basis,dtype=dtype,
                    check_herm=hermitian,check_symm=False,check_pcon=False).toarray())
        hops=[]
        for i in range(L):
            j=(i+1)%L;hops.extend([[-s['t'],i,j],[-s['t'],j,i]])
        H=op([['+-|',hops],['|+-',hops],['n|n',[[s['U'],i,i] for i in range(L)]]])
        nu=[op([['n|',[[1.,j]]]]) for j in range(L)];nd=[op([['|n',[[1.,j]]]]) for j in range(L)]
        up=np.rint(np.diag(sum(nu))).astype(int);down=np.rint(np.diag(sum(nd))).astype(int)
        coords=[(sum(int(round(nu[j][k,k]))<<j for j in range(L)),sum(int(round(nd[j][k,k]))<<j for j in range(L))) for k in range(basis.Ns)]
        return basis,op,H,nu,nd,up,down,coords
    s=cfg['X09'];basis,op,H,nu,nd,up,down,coords=model(s)
    ids=np.flatnonzero((up==s['Nup'])&(down==s['Ndown']));HH=H[np.ix_(ids,ids)];e,v=eigh(HH)
    ground=e-e[0]<=cfg['gates']['ground_cluster_abs'];P=v[:,ground]@v[:,ground].conj().T
    X09={'coordinates':[list(coords[i]) for i in ids], 'Hamiltonian':HH.tolist(),
         'ground_projector':P.real.tolist(),'ground_dimension':int(sum(ground)),
         'ground_energy':float(e[0]),'gap':float(e[~ground][0]-e[0]),
         'eigenvalues':e.tolist(),'eigen_residual':float(np.linalg.norm(HH@v-v*e)),
         'projector_idempotency':float(np.linalg.norm(P@P-P))}
    s=cfg['X12'];L=s['L'];basis,op,H,nu,nd,up,down,coords=model(s)
    ids=np.flatnonzero(up+down==s['N_total']);hc=H[np.ix_(ids,ids)];e,v=eigh(hc)
    keep=e-e[0]<=cfg['gates']['ground_cluster_abs'];V=np.zeros((basis.Ns,int(sum(keep))))
    V[ids,:]=v[:,keep];mixed=V@V.T/sum(keep)
    pureids=np.flatnonzero((up==s['pure_Nup'])&(down==s['pure_Ndown']));ep,vp=eigh(H[np.ix_(pureids,pureids)])
    pu=np.zeros(basis.Ns);pu[pureids]=vp[:,0]
    Sminus=np.zeros_like(H);local_minus=[]
    for j in range(L):
        cu=op([['-|',[[1.,j]]]],False);cd=op([['|-',[[1.,j]]]],False)
        local_minus.append(cd.T@cu);Sminus+=local_minus[-1]
    pd=Sminus@pu;pd/=np.linalg.norm(pd)
    theta=np.pi/s['rotation_theta_pi_denominator'];phi=np.pi/s['rotation_phi_pi_denominator']
    U=np.array([[np.cos(theta),-np.exp(-1j*phi)*np.sin(theta)],
                [np.exp(1j*phi)*np.sin(theta),np.cos(theta)]])
    rotated=V@U;mixrot=rotated@rotated.conj().T/2
    pr=np.cos(theta)*pu+np.exp(1j*phi)*np.sin(theta)*pd
    sz=[(nu[j]-nd[j])/2 for j in range(L)]
    sx=[(m+m.T)/2 for m in local_minus];sy=[(m.T-m)/(2j) for m in local_minus]
    def real(x):
        x=complex(x);assert abs(x.imag)<1e-10
        return float(x.real)
    def measure(label,rho,beta=-1.,epsilon=-1.):
        rdms=[]
        for nA in s['subsystem_sizes']:
            if nA==L:reduced=rho
            else:
                reduced=np.asarray(basis.partial_trace(rho,sub_sys_A=(list(range(nA)),list(range(nA))),return_rdm='A',enforce_pure=False))
                if reduced.shape==(1,4**nA,4**nA):reduced=reduced[0]
                assert reduced.shape==(4**nA,4**nA)
            spectrum=np.linalg.eigvalsh((reduced+reduced.conj().T)/2)[::-1]
            purity=real(np.trace(reduced@reduced))
            rdms.append({'nA':nA,'spectrum':spectrum.tolist(),'trace':real(np.trace(reduced)),
                         'purity':purity,'S2':-math.log(purity)})
        means=[real(np.trace(rho*z)) for z in sz]
        spin_xyz=[[real(np.trace(rho@x)),real(np.trace(rho@y)),real(np.trace(rho@z))] for x,y,z in zip(sx,sy,sz)]
        zz=real(np.trace(rho@sz[0]@sz[1]))
        return {'label':label,'beta':beta,'excited_weight':epsilon,'energy':real(np.trace(rho@H)),
                'trace':real(np.trace(rho)),'Nmean':real(np.trace(rho@sum(nu+nd))),
                'spin_z_means':means,'spin_xyz_means':spin_xyz,'Czz':zz,'connected_Czz':zz-means[0]*means[1],
                'distance_to_equal_ground':float(np.sum(abs(np.linalg.eigvalsh(rho-mixed)))),
                'rdms':rdms}
    rows=[measure('pure_up',np.outer(pu,pu)),measure('pure_rotated',np.outer(pr,pr.conj())),
          measure('equal_ground_mixture',mixed),measure('rotated_basis_ground_mixture',mixrot)]
    for beta in s['beta']:
        weights=np.exp(-beta*(e-e[0]));weights/=sum(weights)
        rhoc=(v*weights)@v.T;rho=np.zeros_like(H);rho[np.ix_(ids,ids)]=rhoc
        rows.append(measure('canonical_gibbs',rho,beta,float(sum(weights[~keep]))))
    X12={'canonical_dimension':len(ids),'canonical_eigenvalues':e.tolist(),
         'ground_dimension':int(sum(keep)),'gap':float(e[~keep][0]-e[0]),
         'eigen_residual':float(np.linalg.norm(hc@v-v*e)),
         'mixture_rotation_error':float(np.linalg.norm(mixrot-mixed)),
         'spin_lowering_norm':float(np.linalg.norm(Sminus@pu)),'rows':rows}
    return {'X09':X09,'X12':X12}

def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',default=str(Path(__file__).with_name('manifest.toml')));p.add_argument('--output',required=True)
    a=p.parse_args();path=Path(a.manifest);cfg=tomllib.loads(path.read_text());start=time.monotonic()
    versions={k:importlib.metadata.version(k) for k in ('quspin','numpy','scipy')}
    assert versions=={k:cfg[k+'_version'] for k in versions}
    result={'schema':cfg['schema'],'status':'COMPLETED','manifest_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'oracle_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'versions':versions}
    result.update(compute(cfg));result['elapsed_seconds']=time.monotonic()-start
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
if __name__=='__main__':main()
