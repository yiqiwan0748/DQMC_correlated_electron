#!/usr/bin/env python3
"""Independent full-Fock continuum/finite-h references; zero scientific completion."""
import argparse,contextlib,datetime,hashlib,importlib.metadata,io,json,math,os,time,tomllib
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'

def compute(cfg):
    import numpy as np
    from scipy.linalg import eigh,expm
    from quspin.basis import spinful_fermion_basis_1d
    from quspin.operators import hamiltonian
    L,U,t=cfg['L'],cfg['U'],cfg['t'];beta=cfg['beta']
    assert (L,U,t,beta,cfg['mu_centered'],cfg['mu_physical'],cfg['boundary'],cfg['nA'])==(4,4.,1.,1.,0.,2.,'pbc',2)
    basis=spinful_fermion_basis_1d(L);assert basis.Ns==4**L
    def operator(static):
        with contextlib.redirect_stdout(io.StringIO()):
            return np.asarray(hamiltonian(static,[],basis=basis,dtype=np.float64,check_herm=True,check_symm=False,check_pcon=False).toarray())
    hops=[]
    for i in range(L):j=(i+1)%L;hops.extend([[-t,i,j],[-t,j,i]])
    K=operator([['+-|',hops],['|+-',hops]])
    Nu=[operator([['n|',[[1.,j]]]]) for j in range(L)]
    Nd=[operator([['|n',[[1.,j]]]]) for j in range(L)]
    D=[operator([['n|n',[[1.,j,j]]]]) for j in range(L)]
    I=np.eye(basis.Ns);N=sum(Nu)+sum(Nd);H=K+U*sum(D)
    V=U*sum(D)-U/2*N+U*L/4*I
    G=K+V;physical_generator=H-cfg['mu_physical']*N+U*L/4*I
    def real(value):
        value=complex(value);assert math.isfinite(value.real) and abs(value.imag)<1e-10
        return float(value.real)
    def maximum(A):return float(np.max(np.abs(A)))
    def reduced(A):
        cut=[x-1 for x in cfg['cut_sites_one_based']]
        result=np.asarray(basis.partial_trace(A,sub_sys_A=(cut,cut),return_rdm='A',enforce_pure=False,sparse=False))
        if result.shape==(1,4**len(cut),4**len(cut)):result=result[0]
        assert result.shape==(4**len(cut),4**len(cut));return result
    def describe(rho):
        A=reduced(rho);purity=real(np.trace(A@A));assert purity>0
        return {'physical_energy':real(np.trace(H@rho)),'Nmean':real(np.trace(N@rho)),
                'purity':purity,'S2':-math.log(purity),'trace':real(np.trace(rho)),'cut_trace':real(np.trace(A))},A
    # Physical local generator, obtained with the same public linear partial trace.
    # The complement dimension removes the trivial trace of identity on B.
    VA=sum(U*(D[j]-.5*(Nu[j]+Nd[j])+.25*I) for j in range(cfg['nA']))
    VA_local=reduced(VA)/(4**(L-cfg['nA']))
    vals,vec=eigh(G);prob=np.exp(-beta*(vals-vals[0]));prob/=sum(prob)
    rho=(vec*prob)@vec.T
    direct=expm(-beta*(G-vals[0]*I));direct/=np.trace(direct)
    continuum,A=describe(rho)
    diagnostics={'physical_mu_centering_matrix_error':maximum(G-physical_generator),
        'K_hermiticity_error':maximum(K-K.T),'V_diagonal_error':maximum(V-np.diag(np.diag(V))),
        'continuum_eigen_residual':maximum(G@vec-vec*vals),
        'continuum_expm_error':maximum(rho-direct),'continuum_rho_hermiticity_error':maximum(rho-rho.T),
        'continuum_rho_eigenvalue_min':float(np.linalg.eigvalsh(rho).min()),
        'continuum_cut_eigenvalue_min':float(np.linalg.eigvalsh(A).min()),
        'continuum_trace_error':abs(continuum['trace']-1),'continuum_cut_trace_error':abs(continuum['cut_trace']-1),
        'continuum_half_filling_error':abs(continuum['Nmean']-L)}
    rows=[];Vdiag=np.diag(V)
    for h in cfg['reference_dtau']:
        M=round(beta/h);assert abs(M*h-beta)<1e-15
        d=np.exp(-h*Vdiag/2);Dhalf=np.diag(d);Dinv=np.diag(1/d)
        kinetic=expm(-h*K);T=np.diag(d*d)@kinetic;S=Dhalf@kinetic@Dhalf
        w,Q=eigh(S);assert min(w)>0
        p=np.exp(M*(np.log(w)-max(np.log(w))));p/=sum(p)
        rhoS=(Q*p)@Q.T;rhoT=Dhalf@rhoS@Dinv
        power=np.linalg.matrix_power(T,M);power/=np.trace(power)
        asym,AT=describe(rhoT);sym,AS=describe(rhoS);DA=expm(-h*VA_local/2)
        transformed=DA@AS@expm(h*VA_local/2)
        row={'dtau':h,'slices':M,'asymmetric':asym,'symmetric':sym,
            'diagnostics':{'transfer_similarity_error':maximum(T-Dhalf@S@Dinv),
                'transfer_eigen_residual':maximum(S@Q-Q*w),'symmetric_transfer_hermiticity_error':maximum(S-S.T),
                'symmetric_transfer_eigenvalue_min':float(w.min()),
                'direct_power_vs_similarity_rho_error':maximum(power-rhoT),
                'full_similarity_rho_error':maximum(rhoT-Dhalf@rhoS@Dinv),
                'cut_similarity_error':maximum(AT-transformed),'cut_purity_similarity_error':abs(asym['purity']-sym['purity']),
                'cut_S2_similarity_error':abs(asym['S2']-sym['S2']),
                'symmetric_rho_eigenvalue_min':float(np.linalg.eigvalsh(rhoS).min()),
                'symmetric_cut_eigenvalue_min':float(np.linalg.eigvalsh(AS).min()),
                'asymmetric_trace_error':abs(asym['trace']-1),'asymmetric_cut_trace_error':abs(asym['cut_trace']-1),
                'symmetric_trace_error':abs(sym['trace']-1),'symmetric_cut_trace_error':abs(sym['cut_trace']-1),
                'half_filling_error':abs(asym['Nmean']-L)}}
        # Do not symmetrize the asymmetric reduced operator: its similarity-invariant purity is native.
        rows.append(row)
    return {'hilbert_dimension':basis.Ns,'continuum':continuum,'diagnostics':diagnostics,'finite_dtau':rows}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--manifest',default=str(Path(__file__).with_name('manifest.toml')));parser.add_argument('--output',required=True)
    args=parser.parse_args();manifest=Path(args.manifest);cfg=tomllib.loads(manifest.read_text());out=Path(args.output)
    assert not out.exists() and cfg['phase']=='PILOT_ONLY' and not cfg['scientific_completion_allowed']
    versions={name:importlib.metadata.version(name) for name in ('quspin','numpy','scipy')}
    for name,version in versions.items():assert version==cfg[name+'_version'],(name,version)
    start=time.time();data=compute(cfg)
    result={'metadata':{'schema':cfg['schema'],'phase':'REFERENCE_ONLY','scientific_completions':0,'catalogue_status':'NOT_RUN',
        'manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(),'adapter_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'versions':versions,'job_id':os.environ.get('SLURM_JOB_ID','NONE'),'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'elapsed_seconds':time.time()-start,'source_urls':cfg['references'],'attribution':'Independent QuSpin full-Fock Hamiltonian/partial trace plus SciPy matrix functions; no native user ED/DQMC import.'},**data}
    out.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')

if __name__=='__main__':main()
