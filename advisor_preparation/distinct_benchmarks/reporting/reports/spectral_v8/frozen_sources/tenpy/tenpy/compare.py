#!/usr/bin/env python3
"""Frozen comparison and eligibility gates; a partial pilot cannot pass an idea."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import tomllib


def padded_error(a,b):
    n=max(len(a),len(b))
    return max((abs((a[i] if i<len(a) else 0.)-(b[i] if i<len(b) else 0.)) for i in range(n)),default=0.)


def spectrum_metrics(a,b,g):
    aa={(s['N'],s['twoSz']):s for s in a}
    bb={(s['N'],s['twoSz']):s for s in b}
    pe,we,xe=0.,0.,0.
    for key in set(aa)|set(bb):
        sa=aa.get(key,{'probabilities':[],'weight':0.})
        sb=bb.get(key,{'probabilities':[],'weight':0.})
        pa,pb=sa['probabilities'],sb['probabilities']
        pe=max(pe,padded_error(pa,pb));we=max(we,abs(sa['weight']-sb['weight']))
        for i,p in enumerate(pa):
            if p>=g['entanglement_energy_probability_min']:
                q=pb[i] if i<len(pb) else 0.
                xe=max(xe,abs(math.log(p)-math.log(q))) if q>0 else 1e100
    return {'probability_abs':pe,'sector_weight_abs':we,'entanglement_energy_abs':xe}


def string_error(a,b):
    aa={(x['i'],x['j']):x['value'] for x in a}
    bb={(x['i'],x['j']):x['value'] for x in b}
    if set(aa)!=set(bb): raise ValueError('Parity interval sets differ')
    return max((abs(aa[k]-bb[k]) for k in aa),default=0.)


def fourier(n):
    import cmath
    L=len(n);mean=sum(n)/L
    return [sum((n[i]-mean)*cmath.exp(-1j*math.pi*m*(i+1)/(L+1)) for i in range(L))/L for m in range(1,L+1)]


def gates(c,o,t,g,idea):
    reasons=[];ed=o['ed']
    if max(ed['residuals'])>g['ed_residual_max']:reasons.append('ED residual gate')
    if ed['gap']<=g['ed_gap_min']:reasons.append('Same-sector gap unresolved')
    for backend,data in [('TeNPy',o),('ITensorMPS',t)]:
        for arm in ('coarse','refined'):
            if arm not in data['arms']:
                reasons.append(f'{backend}/{arm} missing');continue
            x=data['arms'][arm]
            vals=[x[k] for k in ('energy','reported_energy','variance','norm','N','twoSz')]
            if not all(math.isfinite(v) for v in vals):reasons.append(f'{backend}/{arm} nonfinite');continue
            checks=[(abs(x['norm']-1)<=g['norm_abs'],'norm'),
                    (abs(x['N']-c['Nup']-c['Ndown'])<=g['sector_abs'],'N'),
                    (abs(x['twoSz']-c['Nup']+c['Ndown'])<=g['sector_abs'],'twoSz'),
                    (g['variance_min']<=x['variance']<=g['variance_max'],'variance'),
                    (abs(x['energy']-ed['energy'])<=g['energy_vs_ed_abs'],'ED energy'),
                    (abs(x['energy']-ed['energy'])<g['energy_error_fraction_of_gap']*ed['gap'],'gap/energy separation'),
                    (abs(x['reported_energy']-x['energy'])<=g['energy_between_packages_abs'],'reported/contracted energy')]
            reasons += [f'{backend}/{arm}: {label}' for ok,label in checks if not ok]
            if backend=='ITensorMPS' and idea=='P11' and abs(x['rdm_trace']-1)>g['norm_abs']:
                reasons.append(f'{backend}/{arm}: physical prefix trace')
            if backend=='TeNPy' and idea=='P11':
                rd=x['rdm_diagnostics']
                if abs(rd['trace']-1)>g['norm_abs'] or rd['minimum_eigenvalue']< -g['norm_abs']:
                    reasons.append(f'{backend}/{arm}: physical RDM positivity/trace')
                if rd['native_physical_probability_error']>g['native_rdm_probability_abs']:
                    reasons.append(f'{backend}/{arm}: native charge-gauge mapping unresolved')
    if 'refined' in o['arms'] and 'refined' in t['arms']:
        if abs(o['arms']['refined']['energy']-t['arms']['refined']['energy'])>g['energy_between_packages_abs']:
            reasons.append('Independent energies disagree')
    return reasons


def compare(cfg,oracle,target):
    ids=['P11','P12','P15'];rows=[]
    if oracle.get('status')!='COMPLETED' or target.get('status')!='COMPLETED':
        return [{'id':x,'status':'UNRESOLVED','reason':'Incomplete/pilot/failed execution; no scientific completion counted'} for x in ids]
    for case,c in cfg['fixtures'].items():
        o,t=oracle['cases'][case],target['cases'][case]
        g=cfg['gates']
        for idea in c['catalogue_ids']:
            reasons=gates(c,o,t,g,idea)
            if reasons:
                rows.append({'id':idea,'status':'UNRESOLVED','reasons':reasons});continue
            a,b=o['arms']['refined'],t['arms']['refined'];metrics={};stability=0.
            if idea=='P11':
                metrics=spectrum_metrics(a['sectors'],b['sectors'],g)
                for data in (o,t):
                    m=spectrum_metrics(data['arms']['coarse']['sectors'],data['arms']['refined']['sectors'],g)
                    stability=max(stability,m['probability_abs'],m['sector_weight_abs'])
                passed=all(metrics[k]<=g[k] for k in metrics)
            elif idea=='P12':
                fa,fb=fourier(a['density']),fourier(b['density'])
                metrics={'density_abs':padded_error(a['density'],b['density']),
                         'fourier_amplitude_abs':max(abs(x-y) for x,y in zip(fa,fb)),
                         'TeNPy_peak_grid_index':max(range(len(fa)),key=lambda i:abs(fa[i]))+1,
                         'ITensor_peak_grid_index':max(range(len(fb)),key=lambda i:abs(fb[i]))+1}
                for data in (o,t):stability=max(stability,padded_error(data['arms']['coarse']['density'],data['arms']['refined']['density']))
                passed=metrics['density_abs']<=g['density_abs'] and metrics['fourier_amplitude_abs']<=g['fourier_amplitude_abs']
            else:
                metrics={'parity_string_abs':string_error(a['strings'],b['strings'])}
                for data in (o,t):stability=max(stability,string_error(data['arms']['coarse']['strings'],data['arms']['refined']['strings']))
                passed=metrics['parity_string_abs']<=g['parity_string_abs']
            status='UNRESOLVED' if stability>g['convergence_observable_abs'] else ('PASS' if passed else 'FAIL')
            rows.append({'id':idea,'status':status,'metrics':metrics,'coarse_refined_max_probability_or_observable_change':stability,
                         'scope':'Frozen finite-size fixture only; parameter arms and gates are not extra ideas.'})
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',default=str(Path(__file__).with_name('manifest.toml')))
    for k in ('oracle','target','output'):p.add_argument('--'+k,required=True)
    args=p.parse_args();cfg=tomllib.loads(Path(args.manifest).read_text())
    digest=hashlib.sha256(Path(args.manifest).read_bytes()).hexdigest()
    oracle=json.loads(Path(args.oracle).read_text());target=tomllib.loads(Path(args.target).read_text())
    try:
        assert oracle['manifest_sha256']==target['manifest_sha256']==digest,'Manifest receipts differ'
        assert target['target_source_sha256']==cfg['source_sha256'],'Target source differs'
        rows=compare(cfg,oracle,target)
        result={'schema':cfg['schema'],'manifest_sha256':digest,'rows':rows,
                'distinct_ideas_passed':sum(r['status']=='PASS' for r in rows),
                'catalogue_ids':['P11','P12','P15']}
    except Exception as e:
        result={'schema':cfg['schema'],'manifest_sha256':digest,'distinct_ideas_passed':0,
                'rows':[{'id':i,'status':'UNRESOLVED','reason':str(e)} for i in ('P11','P12','P15')]}
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':main()
