#!/usr/bin/env python3
"""Immutable35 spectral_v8 + complete A15/A16/A17 evidence; no physics/jobs."""
import argparse,copy,hashlib,sys,tarfile,tomllib
from pathlib import Path
import leaf_rollup_common_v1 as C

BASE=C.ROOT
BASE35_SHA='25b6733504af033cdbb6d0a942f86bcc8eacbe33c48451602c1af7b5b1f22c23'
BASE35_MANIFEST_SHA='ad8fdb5278f34ac79ea670f48e6449856cdea9d2af1c914cb565c8d5472ef969'
FREEZE_SHA='713582323c240fa37ae4393544d895be76577a975769621aeba850ec069cbd18'
ARCHIVE_SHA='3ac8111fcbf2b04e8201cd4602da80f55b94b672ff335562942671cb3c06d983'
AUDIT_SHA='043c194820b42b6f5d615918a166c95ae210016b68550fec55ab05129f72c9c5'
AUDITOR_SHA='070737427a86b5984bb4f75b82b4bb3b4e26f700e59db9dad85631b2cc968fc9'
RAW_SHA={'oracle.json':'bb245da6b13404459f145b8dbc6e3faefa61ad4b22046dc54b3c421256b20fdd',
         'target.toml':'dd68eec9b659df9e465d252f673a97e689e1acae1b56bafb91f5de1997ed230d',
         'comparison.json':'0e79a01cd4f09d4b35ecf9d501ac5dbfde62a71ff93751f409d811e5e9b5c82a'}
IDS=('A15','A16','A17')

def keyed(rows,fields,wanted):
    keys=[tuple(r[k] for k in fields) for r in rows]
    C.require(len(keys)==len(set(keys))==len(wanted) and set(keys)==set(wanted),'Incomplete/duplicate grid:'+str(fields))

def completeness(m,o,t):
    C.require(C.finite(o) and C.finite(t),'Nonfinite raw receipt')
    C.require(t['mode']=='full' and t['slurm_job_id']=='6458671','Pilot/wrong full job')
    C.require(t['release']==o['release']==m['release']=='two_orbital_analytic_v1','Release identity')
    C.require(t['julia_version']==m['julia_version'] and t['julia_threads']==t['blas_threads']==1,'Target runtime/threads')
    C.require(o['numpy_version']==m['numpy_version'] and o['scipy_version']==m['scipy_version'],'Oracle runtime')
    C.require(t['source_hashes']==m['source_hashes'],'Actual target source receipt')
    wanted=[(u,*s) for u in m['A15']['interactions'] for s in m['A15']['spin_sectors']]
    keyed(t['A15'],('U','Nup','Ndown'),wanted)
    for row in [o['A15'],*t['A15']]:
        C.require(len(row['energies'])==8 and sum(g['multiplicity'] for g in row['groups'])==8,'Incomplete one-electron spectrum')
        for g in row['groups']:
            C.require(len(g['projector'])==64 and len(g['site_probabilities'])==8 and len(g['orbital_weights'])==2 and len(g['energies'])==g['multiplicity'],'Incomplete degenerate projector')
    keyed(t['A16'],('Nk',),[(n,) for n in m['A16']['meshes']])
    width=[(x,) for x in m['A16']['smoothing_widths']]
    for row in [o['A16'],*t['A16']]:
        C.require(len(row['energy'])==m['A16']['bins'],'Partial DOS energy grid')
        keyed(row['curves'],('sigma',),width)
        C.require(all(len(x['density'])==m['A16']['bins'] for x in row['curves']),'Partial DOS curves')
    C.require(len(o['A16']['bin_masses'])==m['A16']['bins'],'Partial continuum bin-mass oracle')
    wanted=[(u,mu,b) for u in m['A17']['interactions'] for mu in m['A17']['centered_chemical_potentials'] for b in m['A17']['inverse_temperatures']]
    for rows in (o['A17'],t['A17']):
        keyed(rows,('U','mu','beta'),wanted)
        C.require(all(len(r['probabilities'])==16 and len(r['spin_connected'])==9 and len(r['local_spin_squared'])==3 for r in rows),'Partial atomic observables')
    return dict(A15_fixed_spin_U_fixtures=6,A15_eigenstates_each=8,A16_meshes=m['A16']['meshes'],A16_widths=m['A16']['smoothing_widths'],A16_bins=512,A17_full_thermal_fixtures=27,A17_fock_states_each=16,all_grids_complete=True)

def verify_verbatim(ev,source,m):
    provenance=C.read(source/'source_provenance.json')
    archive=ev.add('two_orbital_verbatim_baseline_archive',C.PREP/'validation_source/baseline-source.tar.gz',provenance['archive_sha256'])
    matched={}
    with tarfile.open(archive) as tar:
        for r in provenance['files']:
            member=tar.getmember(r['archive_member']);C.require(member.isfile(),'Nonregular baseline source')
            b=tar.extractfile(member).read();f=source/m['source_root']/C.relative(r['path'])
            C.require(hashlib.sha256(b).hexdigest()==r['sha256']==m['source_hashes'][r['path']] and len(b)==r['bytes'] and b==f.read_bytes(),'Baseline production source differs')
            matched[r['path']]=r['sha256']
    C.require(matched==m['source_hashes'],'Incomplete production closure')
    return matched

def gate_metrics(idea,checks):
    expected={'A15':236,'A16':47,'A17':513}[idea]
    C.require(len(checks)==expected and len({r['name'] for r in checks})==expected,'Incomplete/duplicate frozen gate list')
    rows=[]
    for r in checks:
        C.require(r['passed'] is True and isinstance(r['eligibility'],bool),'Failed/unresolved gate')
        q=copy.deepcopy(r)
        if 'bound' not in q:
            C.require(q['value'] is True,'Nonpassing boolean gate');q['relation']='is_true';q['quantity_kind']='boolean eligibility or exact structural condition'
        else:
            ge=q['name'].endswith(('_positivity','_positive_rho','_nonzero_local_spin_moment')) or q['name']=='narrower_broadening_peak_growth'
            q['relation']='greater_than_or_equal' if ge else 'less_than_or_equal'
            C.require((q['value']>=q['bound']) if ge else (q['value']<=q['bound']),'Stored gate does not match inequality')
            q['quantity_kind']='raw signed minimum or peak growth; not absolute error' if ge else 'frozen upper-bound gate value (error, residual or displacement as named)'
        rows.append(q)
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base35',type=Path,default=BASE/'reporting/reports/spectral_v8')
    p.add_argument('--results',type=Path,default=BASE/'results/two_orbital_analytic_full_6458671')
    p.add_argument('--independent-audit',type=Path,required=True)
    p.add_argument('--independent-replay',type=Path,required=True)
    p.add_argument('--independent-auditor',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();C.require(not a.output.exists(),'Choose a new immutable output');a.output.mkdir(parents=True)
    try:
        base,manifest,closure=C.copy_base(a.base35,a.output,BASE35_SHA,BASE35_MANIFEST_SHA,35,'base35');ev=C.Evidence(a.output,manifest)
        ev.add('two_orbital_reporting_tool',__file__);ev.add('two_orbital_reporting_helpers',C.__file__)
        source,freeze=C.frozen_release(ev,'two_orbital',BASE/'source_releases/two-orbital-analytic-v1-freeze.json',FREEZE_SHA,BASE/'source_releases/two-orbital-analytic-v1-source.tar.gz',ARCHIVE_SHA,'two_orbital_analytic',IDS)
        m=tomllib.loads((source/'manifest.toml').read_text());raw={n:ev.add('two_orbital_'+n.replace('.','_'),a.results/n,h) for n,h in RAW_SHA.items()}
        o=C.read(raw['oracle.json']);t=tomllib.loads(raw['target.toml'].read_text());reported=C.read(raw['comparison.json'])
        C.require(o['manifest_sha256']==t['manifest_sha256']==reported['manifest_sha256']==C.sha(source/'manifest.toml'),'Manifest binding')
        C.require(t['driver_sha256']==C.sha(source/'target.jl') and o['driver_sha256']==C.sha(source/'oracle.py') and reported['comparator_sha256']==C.sha(source/'compare.py'),'Actual executable source receipt')
        complete=completeness(m,o,t);matched=verify_verbatim(ev,source,m)
        transfer=C.read(ev.add('two_orbital_remote_transfer',a.results/'remote_transfer_metadata.json'))
        C.require(transfer['job_id']=='6458671' and '6458671|COMPLETED|0:0|' in transfer['sacct'],'Full successful execution receipt')
        for n in ('target.resources.txt','oracle.resources.txt','receipts.sha256','start_utc.txt','end_utc.txt'):
            h=transfer['files']['results/two_orbital_analytic_full_6458671/'+n]['sha256'];ev.add('two_orbital_'+n.replace('.','_'),a.results/n,h)
        for n,path in raw.items():C.require(transfer['files']['results/two_orbital_analytic_full_6458671/'+n]['sha256']==C.sha(path),'Raw transfer binding')
        ev.add('two_orbital_job_log',BASE/'logs/two_orbital_analytic_full_6458671.log',transfer['files']['logs/two_orbital_analytic_full_6458671.log']['sha256'])
        C.require(C.sha(a.independent_audit)==AUDIT_SHA and C.sha(a.independent_auditor)==AUDITOR_SHA and C.sha(a.independent_replay)==RAW_SHA['comparison.json'],'Pinned independent audit/source/replay differs')
        audit,independent=C.independent(ev,'two_orbital',a.independent_audit,a.independent_auditor,a.independent_replay,IDS,raw,ARCHIVE_SHA,'two_orbital_independent_receipt_audit_v1')
        C.require(audit['scientific_completions']==3 and set(audit['benchmarks'])==set(IDS) and all(r['status']=='PASS' for r in audit['benchmarks'].values()),'Partial independent audit')
        result,agreement=C.replay_comparator(ev,'two_orbital',lambda out:[sys.executable,str((source/'compare.py').resolve()),'--target',str(raw['target.toml'].resolve()),'--oracle',str(raw['oracle.json'].resolve()),'--output',str(out.resolve())],reported,independent)
        rmap={r['id']:r for r in result['benchmarks']};C.require(len(result['benchmarks'])==3 and set(rmap)==set(IDS),'Incomplete comparator IDs')
        C.require(all(r['status']=='PASS' for r in rmap.values()) and result['unsupported']==reported['unsupported'] and result['unsupported'][0]['id']=='A25' and result['unsupported'][0]['status']=='NOT_RUN','No partial promotion or A25 inference')
        rows=copy.deepcopy(base['rows'])
        for r in rows:
            idea=r['id']
            if idea not in IDS:continue
            metric=gate_metrics(idea,rmap[idea]['checks']);scope=dict(method='Actual unchanged two-orbital production APIs with explicit new observable adapters versus independent analytical/one-body/product-state oracle',production_DQMC_sampling_exercised=False,production_DMRG_exercised=False,production_sources=matched,
                fixture={'A15':'2x2 open cluster; all8one-electron levels/projectors in six U/spin fixtures','A16':'Square-lattice one-orbital/site/spin DOS from degenerate two-orbital API, four meshes/two fixed broadenings/512bins','A17':'Two independent onsite Hubbard orbitals; full16-state trace at27(U,mu,beta) fixtures'}[idea])
            r.update(status='PASS',reported_status='PASS',reason='',fixtures={'A15':6,'A16':1,'A17':27}[idea],comparisons=len(metric),comparison_unit='complete frozen scientific and eligibility gates; not additional ideas',metrics=metric,scope=scope,input_role='two_orbital_comparison_json',audit_issues=[])
        C.dump(a.output/'analysis/two_orbital_gate_metrics.json',{i:gate_metrics(i,rmap[i]['checks']) for i in IDS})
        C.finish(a.output,base,ev,rows,set(IDS),38,closure,__file__,
            'A15/A16/A17 use unchanged two-orbital Hamiltonian/band/atomic APIs with new full-spectrum, matched-DOS and thermal-trace adapters. They do not establish interacting DQMC sampling, Hund/interorbital interactions, DMRG or arbitrary-flux A25. Original35rows retain their exact scopes.',
            dict(two_orbital=dict(source_archive_sha256=ARCHIVE_SHA,full_job='6458671',completeness=complete,source_closure=matched,independent_receipt_status=audit['status'],comparator_replay=agreement,no_partial_idea_pass=True)),
            'spectral_v8',BASE35_SHA)
    except Exception as ex:
        C.dump(a.output/'ROLLUP_FAILURE.json',dict(status='REJECTED_NO_NEW_LEDGER',reason=f'{type(ex).__name__}: {ex}',physics_executed=False));raise

if __name__=='__main__':main()
