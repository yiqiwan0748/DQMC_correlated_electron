#!/usr/bin/env python3
"""Immutable38 plus exact mixed-status X09/X12 attempt01; no physics/jobs."""
import argparse,copy,csv,hashlib,json,shutil,sys,tomllib
from datetime import datetime,timezone
from pathlib import Path
import leaf_rollup_common_v1 as C

BASE=C.ROOT
BASE38_SHA='ef0ea9e561858941329c28c9f0a78ce9c789d5bc53971b0f1cd7789f24f433d7'
BASE38_MANIFEST_SHA='9da0a99746a2749d06942c13bbc73953eaac85f99dd8ee79864b644eb5761204'
FREEZE_SHA='47d379b14e8d7cf03ee139355f034ec0060694e76ae37c4e868889c63de369d1'
ARCHIVE_SHA='45298d46a32e3f2d2d2af273a989cd8c10f7cd850e7aea631e8c1e6725dc2dd4'
HELPERS_SHA='f2211fcb1fc59c811ccda68587d2abea855b1cde933d9abea97119b064278a38'
IDS=('X09','X12')
# Exact original execution and independently audited mixed-status receipts.
RAW_SHA={'oracle.json':'c41d18eba752608b6b7a415ec3736fbafcc4bff854c267d89ca28c20d6fbde67','target.toml':'f5b2d115cd127522bebb2fdd232f96ab57bfe3ff3d496c361947a60f34f6bb42','comparison.json':'3a5a2d9994918ea135b5a1397334c0eac53bc52fd8c9c7ea4a6d41e9770e9e6a'}
AUDIT_SHA='1568351e7e61800c05ba06532c5f47a7cb4c77c9fdb97959c70f540acd2470b8'
AUDITOR_SHA='5192375591f34ba50e213968233df4d4b6e2da38dfccf80d2842515481a2254f'
AUDIT_SCHEMA='ground_certification_independent_receipt_audit_v1'

def unique(rows,key,wanted):
    ks=[key(r) for r in rows];C.require(len(ks)==len(set(ks))==len(wanted) and set(ks)==set(wanted),'Incomplete/duplicate certification grid')

def completeness(m,o,t):
    C.require(C.finite(o) and C.finite(t),'Nonfinite raw certification')
    C.require(o['schema']==t['schema']==m['schema'] and o['status']==t['status']=='COMPLETED','Incomplete/schema')
    C.require(o['versions']=={k:m[k+'_version'] for k in ('quspin','numpy','scipy')},'Oracle runtime')
    C.require(all(t[k]==m[k] for k in ('julia_version','itensors_version','itensormps_version')),'Target runtime')
    expected={(u,d) for u in range(16) for d in range(16) if u.bit_count()==d.bit_count()==2}
    for r in (o['X09'],t['X09']):
        xy=[tuple(x) for x in r['coordinates']];C.require(len(xy)==len(set(xy))==36 and set(xy)==expected,'Incomplete fixed-spin Hilbert space')
    for owner,r,keys in [('oracle',o['X09'],('Hamiltonian','ground_projector')),('target',t['X09'],('Hamiltonian','MPO_qn_Hamiltonian','MPO_dense_Hamiltonian'))]:
        for k in keys:C.require(len(r[k])==36 and all(len(a)==36 for a in r[k]),'Incomplete '+owner+' matrix')
    p=m['X09'];unique(t['X09']['rows'],lambda r:r['label'],p['starts']+p['challenge_states'])
    for r in t['X09']['rows']:
        C.require(len(r['coeff_real'])==len(r['coeff_imag'])==36,'Partial MPS coefficients')
        if r['label'] in p['starts']:
            n=p['starts'].index(r['label']);C.require(r['kind']=='native_dmrg' and r['seed']==p['seeds'][n] and len(r['truncation_errors'])==len(p['maxdim_schedule']),'Incomplete seeded native DMRG sweeps')
        else:C.require(r['kind']=='manufactured_excited_mps','Misattributed manufactured diagnostic')
    p=m['X12'];wanted=[(x,-1.) for x in p['prescriptions']]+[('canonical_gibbs',x) for x in p['beta']]
    for r in (o['X12'],t['X12']):
        C.require(r['canonical_dimension']==len(r['canonical_eigenvalues'])==4 and r['ground_dimension']==2,'Incomplete canonical dimer manifold')
        unique(r['rows'],lambda x:(x['label'],x['beta']),wanted)
        for a in r['rows']:
            unique(a['rdms'],lambda x:x['nA'],p['subsystem_sizes'])
            C.require(len(a['spin_z_means'])==len(a['spin_xyz_means'])==2 and all(len(v)==3 for v in a['spin_xyz_means']),'Partial physical spin components')
            for z in a['rdms']:C.require(len(z['spectrum'])==4**z['nA'],'Partial subsystem spectrum')
    return dict(X09_fixed_spin_basis=36,X09_native_DMRG_starts=4,X09_manufactured_excited_challenges=2,X12_canonical_states=4,X12_ground_manifold_dimension=2,X12_prescriptions_including_beta=10,X12_RDM_cuts=[1,2],all_grids_complete=True)

def check(rows,name,value,relation,bound,source):
    C.require(C.finite(value) and C.finite(bound),'Nonfinite detailed metric')
    if relation=='less_than_or_equal':passed=value<=bound
    elif relation=='greater_than_or_equal':passed=value>=bound
    elif relation=='inside_closed_interval':passed=bound[0]<=value<=bound[1]
    elif relation=='equals':passed=type(value) is type(bound) and value==bound
    else:raise ValueError('Unknown inequality')
    rows.append(dict(metric=name,value=value,relation=relation,bound=bound,source=source,passed=passed))

def metrics(idea,m,o,t,r):
    rows=[];g=m['gates'];le='less_than_or_equal';ge='greater_than_or_equal';interval='inside_closed_interval'
    for name,value in r['max_abs_errors'].items():
        limit=g['matrix_abs'] if idea=='X09' and name in ('Hamiltonian','MPO_qn_Hamiltonian','MPO_dense_Hamiltonian','MPO_qn_max_imaginary','MPO_dense_max_imaginary') else g['energy_abs'] if idea=='X09' else g['observable_abs']
        check(rows,name,value,le,limit,'frozen comparator absolute-error field; specific bound retained')
    if idea=='X09':
        a=o[idea];b=t[idea]
        check(rows,'unique_external_ground_dimension',a['ground_dimension'],'equals',1,'QuSpin')
        # Native comparator requires a strict gap>; equality is rejected explicitly.
        C.require(a['gap']>g['ground_gap_min'],'External ground gap unresolved')
        check(rows,'external_ground_gap',a['gap'],ge,g['ground_gap_min'],'QuSpin; frozen comparator additionally enforces strict>')
        for name,value in [('oracle_eigen_residual',a['eigen_residual']),('source_eigen_residual',b['source_eigen_residual']),('oracle_projector_idempotency',a['projector_idempotency'])]:check(rows,name,value,le,g['eigen_residual_max'],'raw eligibility')
        cert={x['label']:x for x in r['candidate_certificates']};unique(r['candidate_certificates'],lambda x:x['label'],m['X09']['starts']+m['X09']['challenge_states'])
        for q in b['rows']:
            name=q['label'];v=cert[name];native=q['kind']=='native_dmrg'
            check(rows,name+':variance',v['variance'],interval,[g['variance_min'],g['variance_max']],'raw MPS variance; signed range, not absolute error')
            check(rows,name+':low_variance_gate',v['low_variance_gate_pass'],'equals',True,'frozen comparator certificate')
            check(rows,name+':external_ground_certificate',v['external_ground_certificate_pass'],'equals',native,'frozen comparator certificate; manufactured states must be rejected')
            check(rows,name+':ground_projector_fidelity',v['external_ground_fidelity'],ge if native else le,g['ground_fidelity_min'] if native else g['challenge_fidelity_max'],'independent complete ground projector')
            check(rows,name+':external_energy_excess',v['external_energy_excess'],interval if native else ge,[-g['energy_abs'],g['energy_abs']] if native else .5*a['gap'],'independent ground energy/gap')
            for field,want in [('norm',1.),('outside_sector_weight',0.),('Nup',2.),('Ndown',2.)]:check(rows,name+':'+field+'_deviation',abs(q[field]-want),le,g['norm_abs'],'raw target eligibility')
            if native:check(rows,name+':reported_energy_deviation',abs(q['reported_energy']-q['energy']),le,g['energy_abs'],'native optimizer output versus physical expectation')
            else:check(rows,name+':source_eigen_residual',q['source_eigen_residual'],le,g['eigen_residual_max'],'manufactured diagnostic source')
    else:
        for owner,data in [('QuSpin',o[idea]),('user_ED_thermal_RDM',t[idea])]:
            for key,limit in [('eigen_residual',g['eigen_residual_max']),('mixture_rotation_error',g['matrix_abs'])]:check(rows,owner+':'+key,data[key],le,limit,'raw eligibility')
            for q in data['rows']:
                for z in q['rdms']:
                    name=f"{owner}:{q['label']}:beta={q['beta']}:nA={z['nA']}"
                    check(rows,name+':trace_deviation',abs(z['trace']-1),le,g['rdm_trace_abs'],'complete physical RDM')
                    check(rows,name+':minimum_eigenvalue',min(z['spectrum']),ge,g['rdm_eigenvalue_min'],'signed positivity bound, not absolute error')
            by={(q['label'],q['beta']):q for q in data['rows']};pure=by['pure_up',-1.];mixed=by['equal_ground_mixture',-1.]
            ps={z['nA']:z for z in pure['rdms']};ms={z['nA']:z for z in mixed['rdms']}
            check(rows,owner+':pure_mixed_entropy_separation',abs(ps[1]['S2']-ms[1]['S2']),ge,g['pure_mixed_entropy_separation_min'],'distinct prescribed states')
            check(rows,owner+':pure_mixed_connected_separation',abs(pure['connected_Czz']-mixed['connected_Czz']),ge,g['pure_mixed_connected_separation_min'],'distinct prescribed states')
            check(rows,owner+':last_beta_trace_norm_distance',by['canonical_gibbs',m['X12']['beta'][-1]]['distance_to_equal_ground'],le,g['last_beta_trace_distance_max'],'full canonical Gibbs versus equal ground mixture')
        check(rows,'user_pure_lanczos_residual',t[idea]['pure_lanczos_residual'],le,g['eigen_residual_max'],'actual user Lanczos')
        check(rows,'user_pure_lanczos_energy_error',abs(t[idea]['pure_lanczos_energy']-3),le,g['energy_abs'],'physical dimer convention')
    return rows

def finish_mixed(out,base,ev,rows,closure,rmap,complete,pins,agreement,audit):
    before={r['id']:r for r in base['rows']};preserved={}
    for r in rows:
        if r['id'] not in IDS:
            C.require(r==before[r['id']],'Prior row changed:'+r['id']);preserved[r['id']]=C.canonical(r)
        else:C.require(before[r['id']]['status']=='NOT_RUN' and r['status']==rmap[r['id']]['status'],'Attempt status altered')
    C.require(sum(r==before[r['id']] and r['status']=='PASS' for r in rows)==38,'Prior38PASS changed')
    tally=C.tally(rows);expected=dict(base['counts']);expected['NOT_RUN']-=2
    for r in rmap.values():expected[r['status']]+=1
    C.require(tally==expected,'Mixed status counts disagree')
    receipts={p.relative_to(out).as_posix():C.sha(p) for p in sorted((out/'analysis').rglob('*')) if p.is_file()}
    scope='X12 compares prescribed pure/coherent/ground-manifold/Gibbs states of the canonical interacting dimer. X09 retains the original eight-sweep native DMRG attempt and its unresolved precision gates; manufactured excited states are explicit diagnostics, not natural optimizer failures. No threshold or prior result is changed.'
    data=dict(generated_utc=datetime.now(timezone.utc).isoformat(),tool_sha256=C.sha(__file__),counts=tally,
        category_counts={cat:{s:sum(r['category']==cat and r['status']==s for r in rows) for s in tally} for cat in sorted({r['category'] for r in rows})},
        counting_rule=base['counting_rule'],scope_limit=scope,base_ledger_sha256=BASE38_SHA,
        base_preservation=dict(label='two_orbital_v9',all_prior_PASS_rows_unchanged=True,all_other98row_objects_unchanged=True,canonical_preserved_row_sha256=preserved,all_base_files_copied=len(closure)),
        tranche_audits=dict(ground_attempt01=dict(source_archive_sha256=ARCHIVE_SHA,full_job='6458708',completeness=complete,source_environment_closure=pins,
            independent_receipt_audit_status=audit['status'],scientific_statuses={i:rmap[i]['status'] for i in IDS},comparator_replay=agreement,failed_gates_retained=True,precision_repair_applied=False)),
        input_manifest=ev.manifest,analysis_receipts=receipts,rows=rows)
    C.dump(out/'input_manifest.json',ev.manifest);C.dump(out/'ledger_100.json',data)
    with (out/'ledger_100.csv').open('w',newline='') as f:
        fields=['id','category','title','status','reported_status','fixtures','comparisons','comparison_unit','reason','metrics'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(r[k],sort_keys=True) if k=='metrics' else r.get(k,'') for k in fields})
    (out/'README.md').write_text('# Ground certification: original mixed-status attempt\n\nExactly 100 distinct catalogue IDs: '+str(tally['PASS'])+' PASS, '+str(tally['FAIL'])+' FAIL, '+str(tally['UNRESOLVED'])+' UNRESOLVED, '+str(tally['NOT_RUN'])+' NOT_RUN. X12 passes; X09 remains unresolved at the original frozen precision. All 38 prior passing rows and the other 98 row objects are unchanged.\n\n'+scope+'\n\nAll original matrices, six candidate states, prescribed-state grids, source/runtime closure, comparator and independent audit are retained. Detailed metrics preserve failed variance, fidelity and energy gates with their actual inequality directions. A later precision attempt requires a separate release and output; this record is never overwritten. No physics or cluster job is run by this reporter.\n')
    print(json.dumps(dict(counts=tally,output=str(out),scientific_statuses={i:rmap[i]['status'] for i in IDS},physics_executed=False),indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base38',type=Path,default=BASE/'reporting/reports/two_orbital_v9')
    p.add_argument('--results',type=Path,default=BASE/'results/ground_certification_6458708')
    p.add_argument('--independent-audit',type=Path,required=True);p.add_argument('--independent-replay',type=Path,required=True);p.add_argument('--independent-auditor',type=Path,required=True)
    p.add_argument('--source-root',type=Path,default=C.PREP/'validated_package/HubbardRenyiDQMC');p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    C.require(not a.output.exists(),'Choose a new immutable output');a.output.mkdir(parents=True)
    try:
        C.require(C.sha(C.__file__)==HELPERS_SHA,'Frozen helper source changed')
        base,manifest,closure=C.copy_base(a.base38,a.output,BASE38_SHA,BASE38_MANIFEST_SHA,38,'base38');ev=C.Evidence(a.output,manifest)
        ev.add('ground_reporting_tool',__file__);ev.add('ground_reporting_helpers',C.__file__)
        source,freeze=C.frozen_release(ev,'ground',BASE/'source_releases/ground-certification-v1-freeze.json',FREEZE_SHA,BASE/'source_releases/ground-certification-v1-source.tar.gz',ARCHIVE_SHA,'ground_certification',IDS)
        m=tomllib.loads((source/'manifest.toml').read_text());raw={n:ev.add('ground_'+n.replace('.','_'),a.results/n,h) for n,h in RAW_SHA.items()}
        o=C.read(raw['oracle.json']);t=tomllib.loads(raw['target.toml'].read_text());reported=C.read(raw['comparison.json'])
        C.require(o['manifest_sha256']==t['manifest_sha256']==reported['manifest_sha256']==C.sha(source/'manifest.toml'),'Manifest binding')
        C.require(o['oracle_source_sha256']==C.sha(source/'oracle_quspin.py') and t['adapter_sha256']==C.sha(source/'target.jl') and reported['comparator_sha256']==C.sha(source/'compare.py'),'Actual adapter/comparator source')
        complete=completeness(m,o,t);user=a.output/'frozen_sources/ground_user';pins={}
        for rel,key in [('analysis/compute_s2_ed.jl','ed_source_sha256'),('analysis/audit_u4_l4_thermal_ed.jl','thermal_source_sha256'),('analysis/compute_s2_itensor_dmrg.jl','dmrg_source_sha256'),('Project.toml','project_sha256'),('Manifest.toml','environment_manifest_sha256')]:
            C.require(t[key]==m[key],'Actual source/runtime receipt');f=ev.add('ground_user_'+key,a.source_root/rel,m[key]);q=user/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f,q);pins[rel]=m[key]
        C.dump(a.output/'analysis/ground_user_source_closure.json',pins)
        C.require(C.sha(a.independent_audit)==AUDIT_SHA and C.sha(a.independent_auditor)==AUDITOR_SHA and C.sha(a.independent_replay)==RAW_SHA['comparison.json'],'Pinned independent audit/source/replay')
        audit,independent=C.independent(ev,'ground',a.independent_audit,a.independent_auditor,a.independent_replay,IDS,raw,ARCHIVE_SHA,AUDIT_SCHEMA)
        for n in ('oracle.resources.txt','target.resources.txt','remote_transfer_metadata.json'):ev.add('ground_'+n.replace('.','_'),a.results/n)
        ev.add('ground_job_log',BASE/'logs/ground_certification_6458708.log');ev.add('ground_submission',BASE/'submissions/ground_certification_6458708.json')
        result,agreement=C.replay_comparator(ev,'ground',lambda out:[sys.executable,str((source/'compare.py').resolve()),'--manifest',str((source/'manifest.toml').resolve()),'--oracle',str(raw['oracle.json'].resolve()),'--target',str(raw['target.toml'].resolve()),'--source-root',str(user.resolve()),'--output',str(out.resolve())],reported,independent)
        rmap={r['id']:r for r in result['rows']};C.require(len(result['rows'])==2 and set(rmap)==set(IDS),'Incomplete comparator IDs')
        C.require(result['distinct_ideas_passed']==sum(r['status']=='PASS' for r in rmap.values())==audit['scientific_completions'],'Scientific completion count mismatch')
        C.require(all(audit['benchmarks'][i]['status']==rmap[i]['status'] for i in IDS),'Independent scientific statuses disagree')
        rows=copy.deepcopy(base['rows']);detailed={}
        for r in rows:
            idea=r['id']
            if idea not in IDS:continue
            detailed[idea]=metrics(idea,m,o,t,rmap[idea]);status=rmap[idea]['status']
            if status=='PASS':C.require(rmap[idea]['eligibility_reasons']==[] and all(g['passed'] for g in detailed[idea]),'PASS with failed detailed gate')
            r.update(status=status,reported_status=status,reason='; '.join(rmap[idea]['eligibility_reasons']),fixtures=1,comparisons=6 if idea=='X09' else 20,comparison_unit='four native optimizer starts plus two manufactured challenges' if idea=='X09' else 'ten distinct prescribed ensembles times two RDM cuts; not additional ideas',metrics=detailed[idea],
                scope=dict(method='Actual user Hubbard MPO/ED/Lanczos/Gibbs/RDM functions with explicit new certificate/ensemble adapters versus independent QuSpin',production_DQMC_sampling_exercised=False,production_DMRG_exercised=idea=='X09',source_environment_closure=pins,scientific_status=status,frozen_eight_sweep_budget=idea=='X09',
                    fixture='L4 U4 half-filled periodic Hubbard: four actual DMRG starts tested against the independent certificate; two deliberately manufactured excited states must be rejected' if idea=='X09' else 'Interacting periodic dimer at canonical total N3: pure/coherent/equal-manifold/rotated-mixture/Gibbs prescriptions compared to their own references',
                    limitation='Manufactured rejection does not demonstrate a natural optimizer failure; no two-orbital optimizer certification inferred' if idea=='X09' else 'Exact small-dimer prescribed-state limit; no unique-ground assumption and no native DQMC beta convergence inferred'),input_role='ground_comparison_json',audit_issues=copy.deepcopy(rmap[idea]['eligibility_reasons']))
        C.dump(a.output/'analysis/ground_detailed_metrics.json',detailed)
        finish_mixed(a.output,base,ev,rows,closure,rmap,complete,pins,agreement,audit)

    except Exception as ex:
        C.dump(a.output/'ROLLUP_FAILURE.json',dict(status='REJECTED_NO_NEW_LEDGER',reason=f'{type(ex).__name__}: {ex}',physics_executed=False));raise

if __name__=='__main__':main()
