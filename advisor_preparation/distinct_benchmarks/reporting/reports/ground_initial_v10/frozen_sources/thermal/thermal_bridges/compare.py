#!/usr/bin/env python3
import argparse
import hashlib
import json
import math
from pathlib import Path
import tomllib


def evaluate(cfg,oracle,target):
    g=cfg['gates'];rows=[]
    def finite_payload(value):
        if isinstance(value,dict):return all(finite_payload(x) for x in value.values())
        if isinstance(value,list):return all(finite_payload(x) for x in value)
        return not isinstance(value,float) or math.isfinite(value)
    if not finite_payload(oracle) or not finite_payload(target):raise ValueError('Nonfinite backend receipt')
    global_reasons=[]
    for backend,data in [('QuSpin',oracle),('user_thermal_ED',target)]:
        if data.get('status')!='COMPLETED':global_reasons.append(f'{backend}: execution incomplete')
        for key in ('maximum_trace_error','maximum_rdm_trace_error','hermiticity_error','N_commutator_error'):
            if data['diagnostics'][key]>g['trace_abs']:global_reasons.append(f'{backend}: {key}')
    def cross(a,b,keys):
        errors={k:max(abs(x[k]-y[k]) for x,y in zip(a,b)) for k in keys}
        return errors
    for idea in ('X02','X03','P07'):
        reasons=list(global_reasons);spec=cfg[idea]
        aa,bb=oracle[idea]['rows'],target[idea]['rows']
        if len(aa)!=len(spec['beta']) or len(bb)!=len(spec['beta']):raise ValueError(f'{idea}: grid length differs')
        if [x['beta'] for x in aa]!=spec['beta'] or [x['beta'] for x in bb]!=spec['beta']:raise ValueError(f'{idea}: beta grid differs')
        metrics={};detail={};observed_mismatch=False
        if idea=='X02':
            for backend,data in [('QuSpin',oracle),('user_thermal_ED',target)]:
                z=data[idea]['ground']
                if z['ground_residual']>g['eigen_residual_max'] or z['global_grand_gap']<=g['ground_gap_min']:
                    reasons.append(f'{backend}: unique ground state not resolved')
                mu=spec['mu_physical']
                if not z['mu_minus']+g['charge_gap_margin']<mu<z['mu_plus']-g['charge_gap_margin']:
                    reasons.append(f'{backend}: physical mu not strictly inside charge gap')
                for key,wanted in [('Nmean',spec['N_ground']),('Nup_mean',spec['Nup_ground']),('Ndown_mean',spec['Ndown_ground'])]:
                    if abs(z[key]-wanted)>g['scalar_abs']:reasons.append(f'{backend}: ground {key}')
                for x in data[idea]['rows']:
                    if x['entropy_bound_valid'] and abs(x['delta_S2'])>x['entropy_bias_bound']+g['bound_roundoff_allowance']:
                        reasons.append(f'{backend}: entropy bound violated at beta={x["beta"]}')
                if not any(x['qualifies'] for x in data[idea]['rows']):reasons.append(f'{backend}: no cooling point meets declared bias target')
            z=target[idea]['ground']
            if z['fixed_Nup_Ndown_ground_residual']>g['eigen_residual_max']:reasons.append('User fixed-spin ground residual')
            if abs(z['fixed_Nup_Ndown_ground_S2']-z['S2'])>g['scalar_abs'] or abs(z['fixed_Nup_Ndown_ground_Ephysical']-z['Ephysical'])>g['scalar_abs']:
                observed_mismatch=True
            ground_keys=('Ephysical','S2','purity','Nmean','Nup_mean','Ndown_mean','global_grand_gap','mu_minus','mu_plus')
            metrics.update({'ground_'+k:abs(oracle[idea]['ground'][k]-target[idea]['ground'][k]) for k in ground_keys})
            keys=('mu_physical','mu_centered','S2','purity','delta_S2','Nmean','weight_N_ground','weight_fixed_Nup_Ndown','ground_weight','excited_weight','entropy_bias_bound')
            metrics.update(cross(aa,bb,keys))
            if [(x['entropy_bound_valid'],x['qualifies']) for x in aa]!=[(x['entropy_bound_valid'],x['qualifies']) for x in bb]:observed_mismatch=True
            detail['certified_beta_values_QuSpin']=[x['beta'] for x in aa if x['qualifies']]
            detail['scope']='L4 fixture only; not a bound for the historical L10 CV calculation.'
        elif idea=='X03':
            for backend,data in [('QuSpin',oracle),('user_thermal_ED',target)]:
                for x in data[idea]['rows']:
                    Nc=spec['N_total']
                    if x['mean_N_error']>spec['mean_N_tolerance'] or abs(x['Nmean_GCE']-Nc)>spec['mean_N_tolerance']:
                        reasons.append(f'{backend}: matched-mean root at beta={x["beta"]}')
                    if abs(x['Nmean_canonical']-Nc)>g['scalar_abs'] or abs(x['Nvariance_canonical'])>g['scalar_abs']:
                        reasons.append(f'{backend}: fixed-total-N projector')
                    if x['canonical_dimension']!=math.comb(2*cfg['L'],Nc) or x['fixed_Nup2_Ndown1_dimension']!=24:
                        reasons.append(f'{backend}: fixed-N confused with fixed-spin sector')
                    if x['compressibility']<=spec['minimum_compressibility'] or x['mu_bracket_width']>1.01*spec['mu_bracket_width']:
                        reasons.append(f'{backend}: chemical-potential root unresolved')
            keys=('mu_physical','mu_centered','Nmean_GCE','Nmean_canonical','Nvariance_canonical','compressibility',
                  'canonical_weight_Nup2_Ndown1','Ecanonical_per_site','EGCE_per_site','delta_energy_per_site',
                  'purity_canonical','purity_GCE','delta_purity','S2_canonical','S2_GCE')
            metrics=cross(aa,bb,keys)
            detail['ensemble_bias_QuSpin']=[{'beta':x['beta'],'delta_energy_per_site':x['delta_energy_per_site'],'delta_purity':x['delta_purity']} for x in aa]
            detail['scope']='Each ensemble is checked against its own oracle; their difference is a result, not an error to force to zero.'
        else:
            for backend,data in [('QuSpin',oracle),('user_thermal_ED',target)]:
                separation=max(abs(x['C_mu']-x['beta2_var_grand']) for x in data[idea]['rows'])
                if separation<spec['minimum_formula_separation']:
                    reasons.append(f'{backend}: fixture does not distinguish the two response formulas')
                for x in data[idea]['rows']:
                    if abs(x['mu_physical']-spec['mu_physical'])>g['scalar_abs']:reasons.append(f'{backend}: physical mu retuned')
                    if abs(x['mu_centered']-(spec['mu_physical']-cfg['U']/2))>g['scalar_abs']:reasons.append(f'{backend}: centered/physical mu conversion')
                    if x['finite_difference_error']>g['finite_difference_abs'] or x['covariance_identity_error']>g['scalar_abs']:
                        reasons.append(f'{backend}: fixed-mu response derivative identity')
            keys=('mu_physical','mu_centered','Ephysical','Nmean','Gphysical','E2','N2','G2','EN','EG',
                  'var_E','var_N','var_G','cov_EN','cov_EG','C_mu','beta2_var_grand','C_mu_finite_difference')
            metrics=cross(aa,bb,keys)
            detail['scope']='Fixed-physical-mu d<Ephysical>/dT; neither fixed-N Cv nor grand-generator variance alone.'
        if not all(math.isfinite(x) for x in metrics.values()):reasons.append('Nonfinite comparison')
        tolerance=g['scalar_abs']
        observed_mismatch |= any(x>tolerance for x in metrics.values())
        status='UNRESOLVED' if reasons else ('FAIL' if observed_mismatch else 'PASS')
        rows.append({'id':idea,'status':status,'max_abs_errors':metrics,'eligibility_reasons':reasons,'details':detail})
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',default=str(Path(__file__).with_name('manifest.toml')))
    for k in ('oracle','target','output'):p.add_argument('--'+k,required=True)
    args=p.parse_args();cfg=tomllib.loads(Path(args.manifest).read_text())
    digest=hashlib.sha256(Path(args.manifest).read_bytes()).hexdigest()
    try:
        oracle=json.loads(Path(args.oracle).read_text());target=tomllib.loads(Path(args.target).read_text())
        assert oracle['manifest_sha256']==target['manifest_sha256']==digest,'Manifest receipt mismatch'
        assert oracle['oracle_source_sha256']==hashlib.sha256(Path(__file__).with_name('oracle_quspin.py').read_bytes()).hexdigest(),'Oracle source receipt mismatch'
        assert target['adapter_sha256']==hashlib.sha256(Path(__file__).with_name('target_thermal.jl').read_bytes()).hexdigest(),'Target adapter receipt mismatch'
        for package in ('quspin','numpy','scipy'):assert oracle['backend_versions'][package]==cfg[package+'_version'],package+' version'
        assert target['julia_version']==cfg['julia_version'],'Julia version'
        for key in ('thermal_source_sha256','ed_source_sha256'):assert target[key]==cfg[key],key
        rows=evaluate(cfg,oracle,target)
    except Exception as e:
        rows=[{'id':x,'status':'UNRESOLVED','reason':str(e)} for x in ('X02','X03','P07')]
    result={'schema':cfg['schema'],'manifest_sha256':digest,'rows':rows,'distinct_ideas_passed':sum(x['status']=='PASS' for x in rows)}
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':main()
