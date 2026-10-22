#!/usr/bin/env python3
"""Independent Lieb-Wu thermodynamic integrals; never imports target data/code."""
import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import time
import tomllib
import warnings
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'


def compute(cfg):
    import numpy as np
    from scipy.integrate import quad,IntegrationWarning
    from scipy.special import expit,j0,j1,roots_legendre
    U,t=cfg['U'],cfg['t'];s=cfg['oracle'];X=s['cutoff_x'];a=U/(2*t)
    assert U>0 and t>0
    def integrand(kind,x,u=U):
        f=expit(-u*x/(2*t))
        ratio=j1(x)/x if np.ndim(x)>0 or x!=0 else .5
        if kind=='energy':return -4*t*j0(x)*ratio*f
        if kind=='gap_integral':return 8*t*ratio*f
        return 2*j0(x)*j1(x)*f*(1-f)
    def adaptive(kind,u=U):
        with warnings.catch_warnings():
            warnings.simplefilter('error',IntegrationWarning)
            return quad(lambda x:integrand(kind,x,u),0,X,epsabs=s['quad_epsabs'],
                        epsrel=s['quad_epsrel'],limit=s['quad_limit'])
    values={}
    for kind in ('energy','gap_integral','double_occupancy'):
        q,err=adaptive(kind);gauss={}
        for order in s['gauss_orders']:
            nodes,weights=roots_legendre(order);x=(nodes+1)*X/2
            gauss[str(order)]=float(X/2*np.dot(weights,integrand(kind,x)))
        factor={'energy':4*t,'gap_integral':8*t,'double_occupancy':2}[kind]
        tail=factor*math.exp(-a*X)/(a*X if kind!='double_occupancy' else a)
        estimate=max(err,*[abs(q-v) for v in gauss.values()],s['error_floor'])+tail
        values[kind]={'adaptive_integral':q,'quad_error_estimate':err,'fixed_gauss':gauss,
                      'absolute_tail_bound':tail,'estimated_absolute_error':estimate}
    values['energy']['value']=values['energy']['adaptive_integral']
    values['double_occupancy']['value']=values['double_occupancy']['adaptive_integral']
    values['charge_gap']=dict(values['gap_integral'])
    values['charge_gap']['value']=U-4*t+values['gap_integral']['adaptive_integral']
    h=s['energy_derivative_h']
    E=lambda u:adaptive('energy',u)[0]
    derivative=(E(U-2*h)-8*E(U-h)+8*E(U+h)-E(U+2*h))/(12*h)
    error=abs(derivative-values['double_occupancy']['value'])
    valid=all(values[k]['estimated_absolute_error']<=s['maximum_estimated_absolute_error']
              for k in ('energy','charge_gap','double_occupancy')) and error<=s['energy_derivative_check_abs']
    return {'oracle_gate':'PASS' if valid else 'UNRESOLVED','values':values,
            'energy_derivative_five_point':derivative,'energy_derivative_error':error,
            'error_scope':'Tail bound rigorous using |J_n|<=1; adaptive/Gauss and roundoff errors are numerical estimates, not interval-arithmetic proof.',
            'formulas':{'energy':'-4t integral J0(x)J1(x)/[x(1+exp(Ux/2t))] dx',
                        'charge_gap':'U-4t+8t integral J1(x)/[x(1+exp(Ux/2t))] dx; full mu+ minus mu-',
                        'double_occupancy':'2 integral J0(x)J1(x) f(x)(1-f(x)) dx; f=expit(-Ux/2t)'},
            'sources':['https://arxiv.org/pdf/cond-mat/0207529']}


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',default=str(Path(__file__).with_name('pilot.toml')));p.add_argument('--output',required=True)
    args=p.parse_args();path=Path(args.manifest);cfg=tomllib.loads(path.read_text());start=time.monotonic()
    versions={k:importlib.metadata.version(k) for k in ('numpy','scipy')}
    for k,v in versions.items():assert v==cfg[k+'_version']
    result={'schema':cfg['schema'],'phase':cfg['phase'],'status':'COMPLETED',
            'manifest_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'oracle_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'versions':versions}
    result.update(compute(cfg));result['elapsed_seconds']=time.monotonic()-start
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)

if __name__=='__main__':main()
