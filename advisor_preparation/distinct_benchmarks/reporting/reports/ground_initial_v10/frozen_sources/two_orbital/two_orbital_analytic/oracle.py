#!/usr/bin/env python3
"""Independent analytic/one-body/product-state oracles. No target imports."""
import argparse, hashlib, json, math, pathlib, platform, time, tomllib
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import ellipkm1

HERE = pathlib.Path(__file__).resolve().parent
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def cluster_one_body(t):
    # Independent Kronecker construction, not the target's bond-enumeration loop.
    # Spatial order(y,x), followed by orbital; all three local axes have dimension2.
    t1,t2,t3,t4=t
    eye=np.eye(2);x=np.array([[0.,1.],[1.,0.]]);j=np.array([[0.,1.],[-1.,0.]])
    ax=np.kron(eye,x);ay=np.kron(x,eye)
    return (-np.kron(ax,np.diag([t1,t2]))-np.kron(ay,np.diag([t2,t1]))
            -t3*np.kron(ax@ay,eye)+t4*np.kron(np.kron(j,j),x))

def groups(values,vectors,tol):
    out=[]; start=0
    while start<len(values):
        stop=start+1
        while stop<len(values) and abs(values[stop]-values[start])<=tol: stop+=1
        p=vectors[:,start:stop]@vectors[:,start:stop].conj().T
        out.append(dict(energies=values[start:stop].tolist(),multiplicity=stop-start,
                        projector=p.real.reshape(-1).tolist(),site_probabilities=np.diag(p).real.tolist(),
                        orbital_weights=[float(np.trace(p[a::2,a::2]).real) for a in range(2)]))
        start=stop
    return out

def a15(m):
    h=cluster_one_body(m['hoppings']); es,vs=np.linalg.eigh(h)
    return dict(energies=es.tolist(),groups=groups(es,vs,m['degeneracy_tolerance']),
                residual=float(np.max(np.abs(h@vs-vs*es))),
                orthogonality=float(np.max(np.abs(vs.conj().T@vs-np.eye(8)))))

def square_cdf(e,opts):
    if e<=-4: return 0.0,0.0
    if e>=4: return 1.0,0.0
    def f(k):
        s=e+2*math.cos(k)
        return 0.0 if s<=-2 else 1.0 if s>=2 else .5+math.asin(s/2)/math.pi
    points=[]
    for z in [(-2-e)/2,(2-e)/2]:
        if -1<z<1: points.append(math.acos(z))
    value,error=quad(f,0,math.pi,points=sorted(points),**opts)
    return value/math.pi,error/math.pi

def a16(m):
    edges=np.linspace(-4,4,m['bins']+1); de=edges[1]-edges[0]
    opts=dict(epsabs=m['quad_epsabs'],epsrel=m['quad_epsrel'],limit=m['quad_limit'])
    def rho(e):
        return float(ellipkm1((e/4)**2))/(2*math.pi**2)
    masses=[]; errors=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        value,error=quad(rho,float(lo),float(hi),**opts)
        masses.append(value);errors.append(error)
    cdfs=[square_cdf(float(e),opts) for e in edges]
    cdfm=np.diff([x[0] for x in cdfs]); masses=np.array(masses)
    curves=[]
    for sigma in m['smoothing_widths']:
        radius=min(m['bins']-1,math.ceil(4*sigma/de))
        kernel=np.exp(-.5*(np.arange(-radius,radius+1)*de/sigma)**2);kernel/=kernel.sum()
        density=np.convolve(masses/de,kernel,mode='same'); density/=density.sum()*de
        curves.append(dict(sigma=sigma,density=density.tolist(),peak=float(density.max())))
    return dict(energy=((edges[:-1]+edges[1:])/2).tolist(),bin_masses=masses.tolist(),
                curves=curves,normalization=float(masses.sum()),
                elliptic_vs_reduced_BZ_bin_error=float(np.max(np.abs(masses-cdfm))),
                maximum_reported_quad_error=float(max(errors+[x[1] for x in cdfs])),
                lower_edge_mass=float(masses[edges[1:]<=-m['edge_window_start']].sum()),
                upper_edge_mass=float(masses[edges[:-1]>=m['edge_window_start']].sum()))

def a17(m):
    # Local order(empty,up,down,double); no fermion helpers from target.
    n=np.diag([0.,1.,1.,2.]);d=np.diag([0.,0.,0.,1.])
    sx=np.zeros((4,4),complex);sx[1,2]=sx[2,1]=.5
    sy=np.zeros((4,4),complex);sy[1,2]=-.5j;sy[2,1]=.5j
    sz=np.diag([0.,.5,-.5,0.]); spins=[sx,sy,sz]
    out=[]
    for u in m['interactions']:
        for mu in m['centered_chemical_potentials']:
            for beta in m['inverse_temperatures']:
                e=np.array([u/4,-u/4-mu,-u/4-mu,u/4-2*mu]);w=np.exp(-beta*(e-e.min()));w/=w.sum()
                mean=lambda op: float(np.trace(np.diag(w)@op).real)
                # Product probabilities in TARGET full-Fock bit order, computed by local occupancies.
                probabilities=[]
                for state in range(16):
                    ux=(state>>0)&1;uy=(state>>1)&1;dx=(state>>2)&1;dy=(state>>3)&1
                    probabilities.append(float(w[ux+2*dx]*w[uy+2*dy]))
                local=[mean(s@s) for s in spins]
                out.append(dict(U=u,mu=mu,beta=beta,probabilities=probabilities,
                                density=2*mean(n),double_occupancy=2*mean(d),
                                density_connected=0.,spin_connected=[0.]*9,
                                local_spin_squared=local,local_density_variance=mean(n@n)-mean(n)**2,
                                S_alm_q0=2*mean(sz@sz),delta_alm=0.))
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    output=pathlib.Path(a.output)
    if output.exists(): raise FileExistsError(output)
    m=tomllib.loads((HERE/'manifest.toml').read_text());start=time.time()
    if np.__version__!=m['numpy_version'] or scipy.__version__!=m['scipy_version']: raise RuntimeError('Unpinned numerical versions')
    result=dict(release=m['release'],manifest_sha256=sha(HERE/'manifest.toml'),driver_sha256=sha(__file__),
                backend='independent_onebody_elliptic_product_state',python_version=platform.python_version(),
                numpy_version=np.__version__,scipy_version=scipy.__version__,A15=a15(m['A15']),A16=a16(m['A16']),A17=a17(m['A17']))
    result['elapsed_seconds']=time.time()-start
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')

if __name__=='__main__': main()
