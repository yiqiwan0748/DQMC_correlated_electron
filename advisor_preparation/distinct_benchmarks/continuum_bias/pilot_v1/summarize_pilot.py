#!/usr/bin/env python3
"""Standard-library pilot receipts/eligibility only; no continuum fit or catalogue PASS."""
import argparse,hashlib,json,math,statistics,tomllib
from pathlib import Path

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def jackknife(numerators,denominators,entropy=False):
    B=len(numerators[0]);assert B>=2 and len(numerators)==len(denominators)
    assert all(len(x)==B and all(math.isfinite(y) for y in x) for x in numerators+denominators)
    totals=[(sum(n),sum(d)) for n,d in zip(numerators,denominators)]
    def value(omit=None):
        ratios=[]
        for n,d,(sn,sd) in zip(numerators,denominators,totals):
            dn=sd if omit is None else sd-d[omit];assert dn>0
            ratio=(sn if omit is None else sn-n[omit])/dn
            assert math.isfinite(ratio)
            if entropy:assert ratio>0
            ratios.append(ratio)
        return -sum(math.log(r) for r in ratios) if entropy else ratios[0]
    full=value();loo=[value(i) for i in range(B)];average=statistics.mean(loo)
    return {'plug_in':full,'bias_corrected':B*full-(B-1)*average,
            'standard_error':math.sqrt((B-1)/B*sum((x-average)**2 for x in loo)),
            'bins':B,'leave_one_out_min':min(loo),'leave_one_out_max':max(loo)}

def coarsen(rows):return [[(x[i]+x[i+1])/2 for i in range(0,len(x),2)] for x in rows]

def finite_tree(value):
    if isinstance(value,float):return math.isfinite(value)
    if isinstance(value,dict):return all(finite_tree(x) for x in value.values())
    if isinstance(value,list):return all(finite_tree(x) for x in value)
    return True

def json_safe(value):
    if isinstance(value,float) and not math.isfinite(value):return str(value)
    if isinstance(value,dict):return {k:json_safe(v) for k,v in value.items()}
    if isinstance(value,list):return [json_safe(x) for x in value]
    return value

def summarize(cfg,manifest,target,oracle_path,target_exit):
    here=manifest.resolve().parent;seeds_path=here/cfg['seed_file'];seeds=tomllib.loads(seeds_path.read_text())
    gates=[];issues=[];receipt_hashes={};bundles=[]
    def check(name,ok,detail=None):
        gates.append({'name':name,'eligible':bool(ok),'detail':detail})
        if not ok:issues.append(name)
    def load(path,kind):
        if not path.exists():check(str(path.name)+' exists',False);return {}
        receipt_hashes[path.name]=sha(path)
        try:return tomllib.loads(path.read_text()) if kind=='toml' else json.loads(path.read_text())
        except Exception as ex:check(str(path.name)+' parse',False,str(ex));return {}
    raw=load(target,'toml');oracle=load(oracle_path,'json');meta=raw.get('metadata',{})
    check('all retained target fields finite',finite_tree(raw))
    check('all retained oracle fields finite',finite_tree(oracle))
    expected_sources=cfg['source_hashes']
    checks={'target completed':raw.get('completion',{}).get('status')=='PILOT_ONLY_COMPLETE',
        'target exit zero':target_exit==0,'target scientific completions zero':meta.get('scientific_completions')==0,
        'target phase':meta.get('phase')=='PILOT_ONLY','target manifest pin':meta.get('manifest_sha256')==sha(manifest),
        'target driver pin':meta.get('adapter_sha256')==sha(here/'target_pilot.jl'),
        'target seed pin':meta.get('seeds_sha256')==sha(seeds_path),'target source pins':meta.get('source_hashes')==expected_sources,
        'target Julia pin':meta.get('julia_version')==cfg['julia_version'],
        'target JDQMC pin':meta.get('jdqmcframework_version')==cfg['jdqmcframework_version'],
        'target single threads':meta.get('julia_threads')==meta.get('blas_threads')==1,
        'target exact sampling configuration':meta.get('sampling')=={k:cfg[k] for k in ('pilot_dtau','bundles_per_dtau','burn_sweeps','measured_sweeps','bins','n_stab','sweeps_per_measurement','lambda_indices','lambda_values','delta_lambda','cut_sites_one_based')},
        'seed namespace':seeds.get('namespace')==cfg['seed_namespace']}
    for name,ok in checks.items():check(name,ok)
    om=oracle.get('metadata',{})
    for name,ok in {'oracle phase':om.get('phase')=='REFERENCE_ONLY','oracle completion zero':om.get('scientific_completions')==0,
        'oracle manifest pin':om.get('manifest_sha256')==sha(manifest),'oracle driver pin':om.get('adapter_sha256')==sha(here/'oracle_quspin.py'),
        'oracle versions':om.get('versions')=={x:cfg[x+'_version'] for x in ('quspin','numpy','scipy')},
        'same compute job':meta.get('job_id')==om.get('job_id') and meta.get('job_id') not in (None,'NONE'),
        'full Fock dimension':oracle.get('hilbert_dimension')==256}.items():check(name,ok)
    continuum_names={'physical_mu_centering_matrix_error','K_hermiticity_error','V_diagonal_error',
        'continuum_eigen_residual','continuum_expm_error','continuum_rho_hermiticity_error',
        'continuum_rho_eigenvalue_min','continuum_cut_eigenvalue_min','continuum_trace_error',
        'continuum_cut_trace_error','continuum_half_filling_error'}
    finite_names={'transfer_similarity_error','transfer_eigen_residual','symmetric_transfer_hermiticity_error',
        'symmetric_transfer_eigenvalue_min','direct_power_vs_similarity_rho_error','full_similarity_rho_error',
        'cut_similarity_error','cut_purity_similarity_error','cut_S2_similarity_error',
        'symmetric_rho_eigenvalue_min','symmetric_cut_eigenvalue_min','asymmetric_trace_error',
        'asymmetric_cut_trace_error','symmetric_trace_error','symmetric_cut_trace_error','half_filling_error'}
    def oracle_diagnostics(label,diag,names):
        check(label+' complete declared diagnostic grid',set(diag)==names)
        for name,value in diag.items():
            valid=isinstance(value,(int,float)) and math.isfinite(value)
            if name.endswith('_min'):valid=valid and value>=cfg['gates']['oracle_positive_eigenvalue_min']
            else:valid=valid and abs(value)<=cfg['gates']['oracle_matrix_abs']
            check(label+' '+name,valid,value)
    oracle_diagnostics('continuum',oracle.get('diagnostics',{}),continuum_names)
    reference=oracle.get('finite_dtau',[]);dt=[x.get('dtau') for x in reference]
    check('all four independent finite-h references',len(dt)==4 and sorted(dt)==sorted(cfg['reference_dtau']) and len(set(dt))==4)
    for r in reference:
        check(f'h={r.get("dtau")} integer slices',r.get('slices')==round(cfg['beta']/r['dtau']))
        oracle_diagnostics(f'h={r.get("dtau")}',r.get('diagnostics',{}),finite_names)
    for name,data in [('continuum',oracle.get('continuum',{}))]+[(f'h={r.get("dtau")}',r.get('asymmetric',{})) for r in reference]:
        for key in ('physical_energy','Nmean','purity','S2','trace','cut_trace'):
            value=data.get(key);check(name+' finite '+key,isinstance(value,(float,int)) and math.isfinite(value),value)
        check(name+' positive purity',isinstance(data.get('purity'),(float,int)) and data.get('purity',0)>0,data.get('purity'))
    energy_rows=raw.get('energy_rows',[]);entropy_rows=raw.get('entropy_node_rows',[])
    ekeys=[(r.get('dtau_index'),r.get('bundle')) for r in energy_rows]
    skeys=[(r.get('dtau_index'),r.get('bundle'),r.get('lambda_index')) for r in entropy_rows]
    expected_e={(s['dtau_index'],s['bundle']) for s in seeds['bundles']}
    expected_s={(i,b,k) for i,b in expected_e for k in cfg['lambda_indices']}
    check('complete unique four-energy grid',len(ekeys)==4 and len(set(ekeys))==4 and set(ekeys)==expected_e)
    check('complete unique sixteen-node grid',len(skeys)==16 and len(set(skeys))==16 and set(skeys)==expected_s)
    energies=dict(zip(ekeys,energy_rows));nodes=dict(zip(skeys,entropy_rows));all_deltas=[];api_seconds=[]
    for s in seeds['bundles']:
        i,b=s['dtau_index'],s['bundle'];label=f'dt{i}_bundle{b}'
        E=energies.get((i,b),{});S=[nodes.get((i,b,k),{}) for k in cfg['lambda_indices']]
        ready=True
        for kind,row,k in [('energy',E,None)]+[('lambda'+str(k),row,k) for k,row in zip(cfg['lambda_indices'],S)]:
            name=label+'_'+kind;side=load(target.parent/(name+'.toml'),'toml')
            check(name+' individual receipt matches aggregate',side==row and bool(row))
            ok=row.get('status')=='COMPLETED';check(name+' completed',ok)
            ready=ready and ok
            if not ok:continue
            check(name+' dtau matches',row.get('dtau')==s['dtau'])
            value=row.get('cumulative_deltaG' if k is None else 'native_cumulative_deltaG')
            good=isinstance(value,(int,float)) and math.isfinite(value) and 0<=value<=cfg['gates']['cumulative_deltaG_max']
            check(name+' native cumulative deltaG',good,value)
            if isinstance(value,(int,float)) and math.isfinite(value):all_deltas.append(value)
            elapsed=row.get('elapsed_seconds');valid=isinstance(elapsed,(int,float)) and math.isfinite(elapsed) and elapsed>0
            check(name+' positive runtime',valid,elapsed)
            if valid:api_seconds.append(elapsed)
            check(name+' acceptance',isinstance(row.get('acceptance'),(int,float)) and 0<=row['acceptance']<=1)
            if k is None:
                check(name+' seed',row.get('seed')==s['energy_seed'])
                ds=[row.get('deltaG_after_burnin')]+row.get('cumulative_deltaG_bins',[])
                valid=len(ds)==cfg['bins']+1 and all(isinstance(x,(int,float)) and math.isfinite(x) and x>=0 for x in ds)
                valid=valid and all(a<=z for a,z in zip(ds,ds[1:])) and ds[-1]==value
                check(name+' cumulative history including burn',valid)
            else:
                check(name+' native seeds',row.get('root_seed')==s['entropy_root_seed'] and [row.get('seed1'),row.get('seed2')]==s['node_replica_seeds'][k])
                check(name+' lambda grid',row.get('lambda_index')==k and row.get('lambda_value')==cfg['lambda_values'][k] and row.get('delta_lambda')==cfg['delta_lambda'])
                check(name+' same spatial cut',row.get('cut_sites_one_based')==cfg['cut_sites_one_based'])
            for key in ('numerator_bins','denominator_bins')+(() if k is not None else ('N_numerator_bins',)):
                values=row.get(key,[]);ok=len(values)==cfg['bins'] and all(isinstance(x,(float,int)) and math.isfinite(x) for x in values)
                check(name+' finite '+key,ok);ready=ready and ok
            phase=row.get('denominator_bins',[])
            check(name+' positive-sign phase bins',len(phase)==cfg['bins'] and all(abs(x-1)<=cfg['gates']['positive_phase_abs'] for x in phase))
        if not ready:continue
        try:
            ns=[r['numerator_bins'] for r in S];ds=[r['denominator_bins'] for r in S]
            ej=jackknife([E['numerator_bins']],[E['denominator_bins']]);sj=jackknife(ns,ds,True)
            ej32=jackknife(coarsen([E['numerator_bins']]),coarsen([E['denominator_bins']]))
            sj32=jackknife(coarsen(ns),coarsen(ds),True)
            nj=jackknife([E['N_numerator_bins']],[E['denominator_bins']])
            ref=next(r['asymmetric'] for r in reference if r['dtau']==s['dtau'])
            bundles.append({'dtau_index':i,'dtau':s['dtau'],'bundle':b,'energy':ej,'S2':sj,'Nmean':nj,'energy_32bin':ej32,'S2_32bin':sj32,
                'descriptive_MC_minus_exact_finite_h':{'energy':ej['bias_corrected']-ref['physical_energy'],'S2':sj['bias_corrected']-ref['S2']}})
            check(label+' paired-ratio and full-functional jackknife finite',True)
        except Exception as ex:check(label+' paired-ratio and full-functional jackknife finite',False,str(ex))
    check('all four bundles analyzed',len(bundles)==4)
    groups=[]
    for h in cfg['pilot_dtau']:
        rows=[r for r in bundles if r['dtau']==h]
        if len(rows)!=2:continue
        result={'dtau':h,'bundles':2,'interpretation':'Descriptive two-bundle mean and standard error only; no confidence/equivalence/continuum-fit decision.'}
        for key in ('energy','S2'):
            values=[r[key]['bias_corrected'] for r in rows]
            result[key]={'mean':statistics.mean(values),'standard_error_independent_bundles':statistics.stdev(values)/math.sqrt(2),
                'minus_continuum':statistics.mean(values)-oracle['continuum']['physical_energy' if key=='energy' else 'S2']}
        groups.append(result)
    seconds=sum(api_seconds);pilot_slices=sum(round(cfg['beta']/h) for h in cfg['pilot_dtau']);full_slices=sum(round(cfg['beta']/h) for h in cfg['reference_dtau'])
    projection=[]
    for measurements in (16384,65536,262144):
        scale=8/cfg['bundles_per_dtau']*full_slices/pilot_slices*(2048+measurements)/(cfg['burn_sweeps']+cfg['measured_sweeps'])
        projection.append({'proposed_measurements':measurements,'api_seconds_scaled':seconds*scale})
    return {'schema':'X01_TIMING_PILOT_SUMMARY_V1','status':'TIMING_PILOT_NUMERICAL_ELIGIBLE' if not issues else 'TIMING_PILOT_UNRESOLVED',
        'catalogue_id':'X01','catalogue_status':'NOT_RUN','scientific_completions':0,'production_decision':'NOT_AUTHORIZED_BY_PILOT',
        'interpretation':'Numerical/source/timing eligibility only. Both observables retained; no production precision, continuum fit, adequacy or scientific equivalence gate evaluated.',
        'manifest_sha256':sha(manifest),'seeds_sha256':sha(seeds_path),'summarizer_sha256':sha(__file__),'raw_receipt_sha256':receipt_hashes,
        'gates':gates,'issues':issues,'bundles':bundles,'descriptive_timepoint_means':groups,
        'maximum_native_cumulative_deltaG':max(all_deltas) if all_deltas else None,
        'timing':{'completed_api_calls':len(api_seconds),'sum_recorded_api_seconds':seconds,'target_main_elapsed_seconds':raw.get('completion',{}).get('elapsed_seconds'),
            'unapproved_production_projections':projection,'projection_limitations':'Work-count estimate only; assumes slice/round proportionality, scales cold overhead, omits oracle/loading/summary and variance/mixing uncertainty. No resource guarantee or execution authorization.'}}

def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',default=str(Path(__file__).with_name('manifest.toml')));p.add_argument('--target',required=True);p.add_argument('--oracle',required=True);p.add_argument('--target-exit-code',type=int,required=True);p.add_argument('--output',required=True)
    a=p.parse_args();manifest=Path(a.manifest);cfg=tomllib.loads(manifest.read_text());assert cfg['phase']=='PILOT_ONLY' and cfg['scientific_completion_allowed'] is False
    out=Path(a.output);assert not out.exists();result=summarize(cfg,manifest,Path(a.target),Path(a.oracle),a.target_exit_code)
    out.write_text(json.dumps(json_safe(result),indent=2,sort_keys=True,allow_nan=False)+'\n');print(result['status']+'; X01 NOT_RUN;0 scientific completions.')

if __name__=='__main__':main()
