#!/usr/bin/env python3
"""Audit pilot quality/resources only: every thermodynamic idea remains UNRESOLVED."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import tomllib


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(x):
    if isinstance(x,dict):return all(finite(v) for v in x.values())
    if isinstance(x,list):return all(finite(v) for v in x)
    return not isinstance(x,float) or math.isfinite(x)


def audit(cfg,oracle,target):
    assert finite(oracle) and finite(target),'Nonfinite receipt'
    assert oracle['status']=='COMPLETED' and target['status']=='PILOT_ONLY_COMPLETE','Incomplete pilot'
    assert oracle['phase']==target['phase']==cfg['phase']=='PILOT_ONLY'
    assert target['scientific_completion_allowed'] is False
    assert oracle['versions']=={x:cfg[x+'_version'] for x in ('numpy','scipy')}
    assert target['julia_version']==cfg['julia_version']
    assert target['ITensors_version']==cfg['itensors_version'] and target['ITensorMPS_version']==cfg['itensormps_version']
    assert target['source_sha256']==cfg['source_sha256']
    assert target['project_sha256']==cfg['project_sha256'] and target['environment_manifest_sha256']==cfg['manifest_sha256']
    expected={(L,sector,arm) for L in cfg['L'] for sector in cfg['sectors'] for arm in cfg['arms']}
    rows=target['rows'];lookup={(x['L'],x['sector'],x['arm']):x for x in rows}
    assert len(rows)==len(lookup) and set(lookup)==expected,'Incomplete/duplicate pilot solve'
    reasons=[];g=cfg['gates'];metrics=[]
    if oracle['oracle_gate']!='PASS':reasons.append('Independent quadrature error/stability gate unresolved')
    for row in rows:
        L,sector,armname=row['L'],row['sector'],row['arm'];arm=cfg['arms'][armname];m=L//2
        Nu=m+(sector=='plus');Nd=m-(sector=='minus')
        name=f'L{L}/{sector}/{armname}'
        assert row['Nup_expected']==Nu and row['Ndown_expected']==Nd,name+': wrong odd/even spin sector'
        sid=cfg['sectors'].index(sector)+1
        assert row['seed']==arm['seed_base']+1000*L+10*sid,name+': seed differs'
        checks=[(abs(row['norm']-1)<=g['norm_abs'],'norm'),
                (abs(row['Nup']-Nu)<=g['sector_abs'] and abs(row['Ndown']-Nd)<=g['sector_abs'],'particle/spin projection'),
                (abs(row['energy']-row['reported_energy'])<=g['reported_contracted_energy_abs'],'reported/contracted energy'),
                (g['variance_min']<=row['variance'] and row['variance_per_site']<=arm['variance_per_site_max'],'energy variance'),
                (0<=row['final_truncerr']<=arm['final_truncerr_max'],'final truncation error')]
        if sector=='half':checks.append((g['double_occupancy_min']<=row['double_occupancy_per_site']<=g['double_occupancy_max'],'half-filled double occupancy'))
        reasons += [name+': '+label for good,label in checks if not good]
        for field in ('site_Nup','site_Ndown','site_double_occupancy'):
            assert len(row[field])==L,name+': incomplete site observables'
        assert abs(sum(row['site_Nup'])-row['Nup'])<=1e-9 and abs(sum(row['site_Ndown'])-row['Ndown'])<=1e-9,name+': site-density normalization'
        assert abs(sum(row['site_double_occupancy'])/L-row['double_occupancy_per_site'])<=1e-10,name+': double occupancy mean'
        assert abs(row['variance_per_site']-row['variance']/L)<=1e-10,name+': variance units'
        assert abs(row['energy_per_site']-row['energy']/L)<=1e-10,name+': energy units'
        assert row['sweep_truncerrors'] and len(row['sweep_truncerrors'])==arm['nsweeps'],name+': missing sweep truncation record'
        assert abs(row['sweep_truncerrors'][-1]-row['final_truncerr'])<=1e-14,name+': truncation receipt'
    for L in cfg['L']:
        arms={}
        for armname in cfg['arms']:
            E={s:lookup[L,s,armname]['energy'] for s in cfg['sectors']}
            D=lookup[L,'half',armname]['double_occupancy_per_site']
            gap=E['plus']+E['minus']-2*E['half'];ph=abs(E['plus']-E['minus']-cfg['U'])
            if ph>g['particle_hole_energy_abs']:reasons.append(f'L{L}/{armname}: particle-hole energy identity')
            arms[armname]={'energy_per_site':E['half']/L,'charge_gap':gap,'double_occupancy':D,'particle_hole_error':ph}
        a,b=arms['coarse'],arms['refined']
        for observable,key in [('energy_per_site','coarse_refined_energy_per_site_abs'),('charge_gap','coarse_refined_charge_gap_abs'),('double_occupancy','coarse_refined_double_occupancy_abs')]:
            if abs(a[observable]-b[observable])>g[key]:reasons.append(f'L{L}: {observable} truncation/seed stability')
        metrics.append({'L':L,'arms':arms,'finite_size_offsets_from_thermodynamic_oracle':{
            k:b[k]-oracle['values'][k if k!='energy_per_site' else 'energy']['value']
            for k in ('energy_per_site','charge_gap','double_occupancy')},
            'scope':'Finite-size offsets are diagnostic, never finite-size equality gates.'})
    return {'pilot_quality':'PILOT_QUALITY_PASS' if not reasons else 'PILOT_QUALITY_UNRESOLVED','quality_reasons':reasons,
            'finite_size_diagnostics':metrics,'solves':len(rows),'sampling_seconds_sum':sum(x['elapsed_seconds'] for x in rows),
            'slowest_solve_seconds':max(x['elapsed_seconds'] for x in rows),'largest_bond_observed':max(x['max_bond'] for x in rows),
            'planning_limit':'This L24-only pilot does not establish L32 or full-size feasibility. Review resource receipts before a separately authorized larger bridge or thermodynamic execution release.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',default=str(Path(__file__).with_name('pilot.toml')))
    for key in ('oracle','target','output'):p.add_argument('--'+key,required=True)
    args=p.parse_args();path=Path(args.manifest);cfg=tomllib.loads(path.read_text());sha=digest(path)
    result={'schema':cfg['schema'],'phase':'PILOT_ONLY','manifest_sha256':sha,'analysis_source_sha256':digest(__file__),
            'distinct_ideas_passed':0,'rows':[{'id':x,'status':'UNRESOLVED','reason':'Pilot only; thermodynamic size/truncation/held-out validation not executed.'} for x in cfg['catalogue_ids']]}
    try:
        oracle=json.loads(Path(args.oracle).read_text());target=tomllib.loads(Path(args.target).read_text())
        assert oracle['manifest_sha256']==target['manifest_sha256']==sha,'Manifest receipt differs'
        assert oracle['oracle_source_sha256']==digest(Path(__file__).with_name('oracle_lieb_wu.py')),'Oracle source differs'
        assert target['adapter_sha256']==digest(Path(__file__).with_name('pilot_target.jl')),'Target adapter differs'
        result.update(audit(cfg,oracle,target))
    except Exception as exc:
        result.update(pilot_quality='PILOT_QUALITY_UNRESOLVED',quality_reasons=[f'{type(exc).__name__}: {exc}'])
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':main()
