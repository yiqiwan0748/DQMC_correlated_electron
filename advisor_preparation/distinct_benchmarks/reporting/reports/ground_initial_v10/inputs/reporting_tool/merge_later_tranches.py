#!/usr/bin/env python3
"""Audit additive tranches without running physics or altering the v4 ledger.

Python 3.11+ stdlib. X06's frozen postprocessor additionally requires SciPy.
Every output directory is new, self-contained and includes this exact tool.
"""
import argparse
import copy
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import tomllib

BASE = Path(__file__).resolve().parents[1]
PREP = BASE.parent
BASE_SHA = 'd32231a48f54a46d5233af6df9d7c0d275adec4f0d0134ee8cf4f2f21d696de0'
RELEASES = {
    'thermal': ('thermal-bridges-v1', 'thermal-bridges-v1', 'thermal_bridges',
                '73e34d1c706f302b79523559bdac0faf7b4c84a62ed9e840f153ad15ac29d826', ['P07', 'X02', 'X03']),
    'tenpy': ('tenpy-tranche', 'tenpy-v1', 'tenpy',
              'cc00d51fc99cf2d265a1fe2c75f59e649a35ff3116c69211511bb83f4987bbed', ['P11', 'P12', 'P15']),
    'coverage': ('dqmc-coverage-pilot-v2', 'dqmc-coverage-pilot-v2', 'dqmc_coverage',
                 '70750c7ef0d9409b63f48cc5678d3c160d1b092e041c28ed58ca749c426dcf1d', ['X06']),
}
STATUSES = {'PASS', 'FAIL', 'UNRESOLVED', 'NOT_RUN'}
PROJECT_SHA = 'ea3b7eacfcf9a36e4652654a5e7b573ef27529edfbad2e272da2bafca1a7f58f'
MANIFEST_SHA = '5cefaf1acdb4fa23a456c2c58f20fb1fcebd54b3cd27bc448cb65b0fcfecb90e'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read(path):
    path = Path(path)
    return json.loads(path.read_text()) if path.suffix == '.json' else tomllib.loads(path.read_text())


def finite(value):
    if isinstance(value, dict):
        return all(finite(x) for x in value.values())
    if isinstance(value, list):
        return all(finite(x) for x in value)
    return not isinstance(value, float) or math.isfinite(value)


def compare_replay(reported, recomputed):
    """Same frozen scientific gates; permit only platform rounding in receipts."""
    differences=[]
    def visit(a, b, path):
        require(type(a) is type(b), f'Replay value type differs: {path}')
        if isinstance(a, dict):
            require(a.keys() == b.keys(), f'Replay keys differ: {path}')
            for key in a:
                visit(a[key], b[key], path+'/'+str(key))
        elif isinstance(a, list):
            require(len(a) == len(b), f'Replay list length differs: {path}')
            for i, (x,y) in enumerate(zip(a,b)):
                visit(x,y,path+'/'+str(i))
        elif isinstance(a, float):
            delta=abs(a-b)
            require(math.isfinite(a) and math.isfinite(b) and
                    delta <= 64*sys.float_info.epsilon*max(1.,abs(a),abs(b)),
                    f'Replay float exceeds machine-roundoff allowance: {path}')
            if delta:
                differences.append({'path':path,'reported':a,'recomputed':b,'absolute_difference':delta})
        else:
            require(a == b, f'Replay value differs: {path}')
    visit(reported,recomputed,'')
    return {'float_allowance':'64*binary64_epsilon*max(1,abs(reported),abs(recomputed)); scientific gates unchanged',
            'nonfloat_values':'exact','roundoff_differences':differences,
            'maximum_absolute_roundoff_difference':max((x['absolute_difference'] for x in differences),default=0.)}


class Evidence:
    def __init__(self, output):
        self.output = output
        self.manifest = {}

    def add(self, role, path, expected=None):
        path = Path(path).resolve()
        require(path.is_file(), f'{role}: missing file {path.name}')
        digest = sha(path)
        if expected:
            require(digest == expected, f'{role}: SHA256 mismatch')
        relative = Path('inputs') / role / path.name
        dest = self.output / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        require(not dest.exists(), f'Duplicate evidence role {role}')
        shutil.copyfile(path, dest)
        try:
            original = path.relative_to(PREP.resolve()).as_posix()
        except ValueError:
            original = f'external/{role}/{path.name}'
        self.manifest[role] = dict(original_path=original, evidence_copy=relative.as_posix(),
                                   sha256=digest, bytes=dest.stat().st_size)
        return dest


def frozen_release(kind, ev):
    stem, release, folder, archive_sha, ids = RELEASES[kind]
    freeze = read(ev.add(kind + '_freeze', BASE / 'source_releases' / (stem + '-freeze.json')))
    archive = ev.add(kind + '_archive', BASE / 'source_releases' / (stem + '-source.tar.gz'), archive_sha)
    require(freeze['release'] == release and freeze['archive_sha256'] == archive_sha, 'Release identity differs')
    require(sorted(freeze.get('distinct_ids', freeze.get('catalogue_ids', []))) == sorted(ids), 'Release ID set differs')
    seen = {}
    with tarfile.open(archive) as tar:
        for member in tar:
            if member.isdir():
                continue
            name = PurePosixPath(member.name)
            require(member.isfile() and not name.is_absolute() and '..' not in name.parts, 'Unsafe archive member')
            require(name.parts[0] == 'distinct_benchmarks', 'Unexpected archive prefix')
            relative = PurePosixPath(*name.parts[1:]).as_posix()
            require(relative in freeze['files'] and relative not in seen, 'Unexpected/duplicate archive file')
            data = tar.extractfile(member).read()
            digest = hashlib.sha256(data).hexdigest()
            require(digest == freeze['files'][relative], f'Frozen file mismatch: {relative}')
            seen[relative] = digest
            dest = ev.output / 'frozen_sources' / kind / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
    require(seen == freeze['files'], 'Archive closure incomplete')
    return ev.output / 'frozen_sources' / kind / folder, freeze


def source_pair(ev, kind, source_root, files):
    for relative, digest in files.items():
        ev.add(kind + '_user_source_' + relative.replace('/', '_'), source_root / relative, digest)


def run_analysis(ev, kind, command):
    proc = subprocess.run(command, cwd=ev.output, capture_output=True, text=True, timeout=600)
    # Logs are analysis evidence, not console output from a physics calculation.
    log = ev.output / 'analysis' / (kind + '_process.json')
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(json.dumps(dict(returncode=proc.returncode, stdout=proc.stdout,
                                  stderr=proc.stderr), indent=2))
    require(proc.returncode in (0, 1), f'{kind}: postprocessor failed (see analysis process receipt)')
    return proc


def check_thermal(cfg, oracle, target):
    require(finite(oracle) and finite(target), 'Nonfinite thermal data')
    for data in (oracle, target):
        require(data['schema'] == cfg['schema'] and data['status'] == 'COMPLETED', 'Thermal execution incomplete')
        for idea in ('P07', 'X02', 'X03'):
            require([x['beta'] for x in data[idea]['rows']] == cfg[idea]['beta'], idea + ': incomplete grid')
        ground = data['X02']['ground']
        for row in data['X02']['rows']:
            eps = row['excited_weight']
            require(0 <= eps <= 1 and abs(eps + row['ground_weight'] - 1) < 1e-10, 'Excited-state weights')
            valid = cfg['X02']['purity_lipschitz_factor'] * eps < ground['purity']
            require(row['entropy_bound_valid'] == valid, 'Entropy bound validity flag')
            bound = -math.log1p(-4 * eps / ground['purity']) if valid else -1.0
            require(abs(row['entropy_bias_bound'] - bound) <= 2e-10, 'Entropy bound formula')
            require(abs(row['delta_S2'] - (row['S2'] - ground['S2'])) < 1e-9, 'Entropy bias identity')
            qualify = (valid and bound <= cfg['X02']['entropy_bias_target'] and
                       max(0., 1 - row['weight_N_ground']) <= cfg['X02']['particle_sector_leakage_target'] and
                       abs(row['delta_S2']) <= cfg['X02']['entropy_bias_target'])
            require(row['qualifies'] == qualify, 'Cooling qualification flag')
            require(abs(row['purity'] - math.exp(-row['S2'])) < 1e-9, 'Purity/entropy identity')
            require(abs(row['purity']-ground['purity']) <= 4*eps+2e-10, 'Purity Lipschitz bound')
        for row in data['P07']['rows']:
            b, mu = row['beta'], row['mu_physical']
            require(abs(row['C_mu'] - b*b*(row['var_E']-mu*row['cov_EN'])) < 1e-8, 'Fixed-mu covariance identity')
            require(abs(row['finite_difference_error'] - abs(row['C_mu']-row['C_mu_finite_difference'])) < 1e-10,
                    'Finite-difference error receipt')


def check_tenpy(cfg, oracle, target):
    require(finite(oracle) and finite(target), 'Nonfinite TeNPy/ITensor data')
    for data in (oracle, target):
        require(data['schema'] == cfg['schema'] and data['status'] == 'COMPLETED' and data.get('pilot') is False,
                'Pilot/incomplete TeNPy execution cannot pass')
        require(set(data['cases']) == set(cfg['fixtures']), 'Missing/extra TeNPy fixture')
        for fixture, spec in cfg['fixtures'].items():
            case = data['cases'][fixture]
            require(set(case['arms']) == set(cfg['arms']), 'Missing/extra convergence arm')
            for arm in case['arms'].values():
                require(len(arm['density']) == spec['L'], 'Incomplete site-density vector')
                if 'P15' in spec['catalogue_ids']:
                    wanted = {(i, j) for i in range(1, spec['L']+1) for j in range(i, spec['L']+1)}
                    actual = [(x['i'], x['j']) for x in arm['strings']]
                    require(len(actual) == len(wanted) and set(actual) == wanted, 'Incomplete/duplicate parity intervals')
                if 'P11' in spec['catalogue_ids']:
                    sectors = arm['sectors']
                    expected = {(up+dn, up-dn) for up in range(spec['cut']+1) for dn in range(spec['cut']+1)}
                    require(len(sectors) == len(expected) and {(s['N'],s['twoSz']) for s in sectors} == expected,
                            'Incomplete/duplicate physical charge sectors')
                    require(abs(sum(s['weight'] for s in sectors)-1) <= 1e-9, 'Sector weights not normalized')
                    require(sum(len(s['probabilities']) for s in sectors) > 0, 'Empty Schmidt spectrum')
                    for sec in sectors:
                        p = sec['probabilities']
                        require(p == sorted(p, reverse=True) and all(x > cfg['schmidt_probability_floor'] for x in p),
                                'Invalid retained Schmidt probabilities')
                        require(sec['weight'] >= -1e-12 and abs(sum(p)-sec['weight']) <= 1e-9, 'Sector probability sum')
                require(abs(sum(arm['density'])-arm['N']) <= 1e-8, 'Density/particle-number identity')
    for fixture, spec in cfg['fixtures'].items():
        ed = oracle['cases'][fixture]['ed']
        require(len(ed['residuals']) == 2 and ed['dimension'] == math.comb(spec['L'],spec['Nup'])*math.comb(spec['L'],spec['Ndown']),
                'Independent ground/gap ED dimension or residual count')


def compare_tranche(kind, directory, ev, source_root):
    own, freeze = frozen_release(kind, ev)
    cfg = read(own / 'manifest.toml')
    oracle_path = ev.add(kind + '_oracle', directory / 'oracle.json')
    target_path = ev.add(kind + '_target', directory / 'target.toml')
    reported_path = (ev.add(kind+'_reported_comparison', directory/'comparison.json')
                     if (directory/'comparison.json').is_file() else None)
    oracle, target = read(oracle_path), read(target_path)
    digest = sha(own / 'manifest.toml')
    require(oracle['manifest_sha256'] == target['manifest_sha256'] == digest, 'Manifest digest differs')
    if kind == 'thermal':
        sources = {'analysis/compute_s2_ed.jl': cfg['ed_source_sha256'],
                   'analysis/audit_u4_l4_thermal_ed.jl': cfg['thermal_source_sha256']}
        require(oracle['oracle_source_sha256'] == sha(own/'oracle_quspin.py'), 'Oracle source differs')
        require(target['adapter_sha256'] == sha(own/'target_thermal.jl'), 'Thermal adapter differs')
        require(all(oracle['backend_versions'][p] == cfg[p+'_version'] for p in ('quspin','numpy','scipy')), 'Oracle versions differ')
        require(target['julia_version'] == cfg['julia_version'], 'Julia version differs')
        require(all(target[k] == cfg[k] for k in ('ed_source_sha256','thermal_source_sha256')), 'User thermal source differs')
        check_thermal(cfg, oracle, target)
    else:
        sources = {'analysis/compute_s2_itensor_dmrg.jl': cfg['source_sha256'],
                   'Project.toml':PROJECT_SHA, 'Manifest.toml':MANIFEST_SHA}
        require(oracle['oracle_source_sha256'] == sha(own/'oracle_tenpy.py'), 'TeNPy adapter differs')
        require(target['adapter_sha256'] == sha(own/'target_itensor.jl'), 'ITensor adapter differs')
        require(target['target_source_sha256'] == cfg['source_sha256'], 'User DMRG source differs')
        require(oracle['versions'] == {'physics-tenpy':cfg['tenpy_version'],'numpy':'2.4.6','scipy':'1.17.1'}, 'TeNPy versions differ')
        require(target['julia_version'] == '1.12.1', 'Julia version differs')
        require(sha(source_root/'Manifest.toml') == MANIFEST_SHA, 'Pinned Julia dependency manifest differs')
        pinned = read(source_root/'Manifest.toml')['deps']
        require(target['library_versions'] == {p:pinned[p][0]['version'] for p in ('ITensorMPS','ITensors')},
                'Runtime ITensor library versions differ from pinned manifest')
        check_tenpy(cfg, oracle, target)
    source_pair(ev, kind, source_root, sources)
    out = ev.output/'analysis'/(kind+'_recomputed.json')
    out.parent.mkdir(parents=True, exist_ok=True)
    proc = run_analysis(ev, kind, [sys.executable, str(own/'compare.py'), '--manifest', str(own/'manifest.toml'),
                 '--oracle', str(oracle_path), '--target', str(target_path), '--output', str(out)])
    require(out.is_file(), f'{kind}: no analysis receipt; ' + proc.stderr.strip().split('\n')[-1])
    result = read(out)
    require(sorted(x['id'] for x in result['rows']) == sorted(RELEASES[kind][4]), 'Comparator ID set differs')
    if reported_path:
        reported = read(reported_path)
        result['reporting_replay']=compare_replay(reported,result)
    return result, freeze


def coverage_tranche(paths, mode, ev, source_root, oracle_path, python):
    own, freeze = frozen_release('coverage', ev)
    cfg = read(own/'manifest.toml')
    source_pair(ev, 'coverage', source_root, read(own/'source_hashes.toml')['files'])
    ev.add('coverage_external_oracle', oracle_path, cfg['oracle_sha256'])
    require(read(oracle_path)['values'][cfg['oracle_case']] == cfg['oracle_reference'], 'Independent oracle value differs')
    copied = [ev.add(f'coverage_raw_{i:03d}', path) for i,path in enumerate(paths,1)]
    # Check the completion sentinel as well as the frozen auditor's raw-bin and range audit.
    for path in copied:
        d = read(path)
        require(d['metadata']['mode'] == mode, 'Pilot/final mode differs')
        require(d['completion']['status'] == ('PILOT_ONLY' if mode == 'pilot' else 'SHARD_COMPLETE'), 'Incomplete coverage shard')
    out = ev.output/'analysis'/'coverage_recomputed.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    proc = run_analysis(ev, 'coverage', [python, str(own/'summarize_coverage.py'), '--mode', mode,
                 '--inputs', *[str(p.relative_to(ev.output)) for p in copied], '--output', str(out.relative_to(ev.output))])
    require(out.is_file(), 'coverage: no analysis receipt; ' + proc.stderr.strip().split('\n')[-1])
    result = read(out)
    require(result['analysis_source_sha256'] == sha(own/'summarize_coverage.py'), 'Coverage analysis source differs')
    require(result['manifest_sha256'] == sha(own/'manifest.toml') and result['mode'] == mode, 'Coverage analysis manifest/mode differs')
    if mode == 'pilot':
        result['reported_status'] = result['status']
        result['status'] = 'UNRESOLVED'
        result['reason'] = 'Pilot only; no X06 scientific completion counted'
    elif result['status'] in ('PASS','FAIL'):
        require(result.get('raw_audit') == 'PASS' and result['total_replicates'] == cfg['final_replicates'] == 400,
                'Coverage final400 audit incomplete')
        require([r['id'] for r in result['replicates']] == list(range(1,401)), 'Coverage replicate ID set differs')
    return {'rows':[result]}, freeze


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base-ledger', type=Path, default=BASE/'reporting/reports/tranche01_repaired_v4/ledger_100.json')
    p.add_argument('--thermal-dir', type=Path)
    p.add_argument('--tenpy-dir', type=Path)
    group=p.add_mutually_exclusive_group()
    group.add_argument('--coverage-final', type=Path, nargs='+')
    group.add_argument('--coverage-pilot', type=Path, nargs='+')
    p.add_argument('--coverage-python', default=sys.executable, help='Python executable with SciPy for the frozen X06 auditor')
    p.add_argument('--coverage-oracle', type=Path, default=PREP/'isaac/results/quspin_attempt02/oracle_values.toml')
    p.add_argument('--renyi-source-root', type=Path, default=PREP/'validated_package/HubbardRenyiDQMC')
    p.add_argument('--history', type=Path, nargs='*', default=[], help='Exact predecessor receipts; never promoted to PASS')
    p.add_argument('--output-dir', type=Path, required=True)
    args=p.parse_args()
    output=args.output_dir.resolve()
    require(not output.exists(), 'Output revision already exists; choose a new directory')
    require(sha(args.base_ledger) == BASE_SHA, 'Base ledger is not the independently audited v4 ledger')
    base=read(args.base_ledger)
    require(base['counts'] == {'PASS':25,'FAIL':0,'UNRESOLVED':0,'NOT_RUN':75}, 'Unexpected v4 counts')
    require(len(base['rows']) == len({r['id'] for r in base['rows']}) == 100, 'Base must have exactly100 unique IDs')
    output.mkdir(parents=True)
    ev=Evidence(output)
    ev.add('reporting_tool', Path(__file__))
    ev.add('base_ledger', args.base_ledger, BASE_SHA)
    # Preserve exact previous evidence, including scientific failures and source amendments.
    for role, item in base['input_manifest'].items():
        copied=ev.add('base_'+role, args.base_ledger.parent/item['evidence_copy'], item['sha256'])
        require(copied.stat().st_size == item['bytes'], 'Base evidence byte count differs')
    for i, path in enumerate(args.history,1):
        ev.add(f'predecessor_{i:03d}', path)
    rows=copy.deepcopy(base['rows']); byid={r['id']:r for r in rows}
    audits={}
    supplied={'thermal':args.thermal_dir, 'tenpy':args.tenpy_dir,
              'coverage':args.coverage_final or args.coverage_pilot}
    for kind, source in supplied.items():
        if not source:
            continue
        ids=RELEASES[kind][4]
        try:
            if kind == 'coverage':
                result, freeze=coverage_tranche(source, 'final' if args.coverage_final else 'pilot', ev,
                   args.renyi_source_root, args.coverage_oracle, args.coverage_python)
            else:
                result, freeze=compare_tranche(kind, source.resolve(), ev, args.renyi_source_root)
            audits[kind]={'status':'AUDITED','release':freeze['release'], 'archive_sha256':freeze['archive_sha256']}
            if 'reporting_replay' in result:
                audits[kind]['reporting_replay']=result['reporting_replay']
            if kind == 'coverage' and result['rows'][0].get('raw_audit') != 'PASS':
                audits[kind]['status']='INCOMPLETE_RAW_AUDIT'
                audits[kind]['reason']=result['rows'][0].get('reason','Frozen auditor did not establish a complete raw-bin audit')
        except Exception as exc:
            reason=f'{type(exc).__name__}: {exc}'
            audits[kind]={'status':'UNRESOLVED','reason':reason}
            result={'rows':[{'id':i,'status':'UNRESOLVED','reason':reason} for i in ids]}
        for new in result['rows']:
            row=byid[new['id']]
            require(row['status'] == 'NOT_RUN', 'Additive tool cannot silently replace an already executed base idea')
            row['history']=[{'status':row['status'],'reason':row.get('reason',''),'source':'immutable_v4'}]
            require(new['status'] in STATUSES, 'Unknown status')
            row['status']=new['status'];row['reported_status']=new.get('reported_status',new['status'])
            row['reason']=new.get('reason','; '.join(new.get('eligibility_reasons',new.get('reasons',[]))))
            row['input_role']=kind
            row['audit_issues']=[] if audits[kind]['status']=='AUDITED' else [audits[kind]['reason']]
            row['later_tranche_result']=new
            row['metrics']=[{'metric':k,('maximum_error' if kind=='thermal' or k.endswith('_abs') else 'value'):v}
                            for k,v in new.get('max_abs_errors',new.get('metrics',{})).items()]
            if kind=='coverage':
                row['metrics']=[{'metric':k,'value':new[k]} for k in ('observed_coverage','binomial_pvalue','coverage_clopper_pearson99',
                    'mean_bias','mean_bias_interval99') if k in new]
            row['scope']={'method':{'thermal':'Actual user thermal ED/ground ED plus ensemble and response adapters versus QuSpin',
                    'tenpy':'Actual user Hubbard MPO and ITensorMPS DMRG versus independent TeNPy DMRG; new observables adapter',
                    'coverage':'Actual corrected user DQMC sampler and eight-chain estimator versus finite-dtau QuSpin oracle'}[kind],
                'execution_evidence_audited':audits[kind]['status']=='AUDITED',
                'production_DQMC_sampling_exercised':kind=='coverage' and audits[kind]['status']=='AUDITED',
                'production_DMRG_exercised':kind=='tenpy' and audits[kind]['status']=='AUDITED',
                'parameter_points_and_convergence_arms_are_extra_ideas':False}
            row['fixtures']=None;row['comparisons']=None
            row['executed_fixture_scope']={'thermal':'Frozen L4 U4 PBC; exact beta grids and ensembles in copied manifest',
                'tenpy':'Full requires L8 OBC doped/half fixtures and both coarse/refined arms; pilot cannot PASS',
                'coverage':'Final requires exactly400 independent eight-chain replicates; pilot cannot PASS'}[kind]
    counts={s:sum(r['status']==s for r in rows) for s in ('PASS','FAIL','UNRESOLVED','NOT_RUN')}
    categories={c:{s:sum(r['category']==c and r['status']==s for r in rows) for s in counts} for c in sorted({r['category'] for r in rows})}
    ledger={'generated_utc':datetime.now(timezone.utc).isoformat(), 'tool_sha256':sha(__file__),
        'base_ledger_sha256':BASE_SHA,'counts':counts,'category_counts':categories,
        'counting_rule':base['counting_rule'],'scope_limit':base['scope_limit'],
        'base_history_and_source_amendments':'Exact immutable base ledger and every base evidence file copied under inputs; no replacement or deletion',
        'later_tranche_audits':audits,'input_manifest':ev.manifest,
        'analysis_receipts':{p.relative_to(output).as_posix():sha(p) for p in sorted((output/'analysis').glob('*')) if p.is_file()},
        'rows':rows}
    (output/'ledger_100.json').write_text(json.dumps(ledger,indent=2,allow_nan=False)+'\n')
    (output/'input_manifest.json').write_text(json.dumps(ev.manifest,indent=2)+'\n')
    with (output/'ledger_100.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['id','category','title','status','reason','input_role','metrics'])
        w.writeheader()
        for row in rows:
            w.writerow({k:json.dumps(row[k],sort_keys=True) if k=='metrics' else row.get(k,'') for k in w.fieldnames})
    text=['# Read-only later-tranche audit','',', '.join(f'{v} {k}' for k,v in counts.items())+'.',
        '', 'Exactly100 catalogue IDs; fixtures, seeds, replicas, convergence arms and matrix elements do not add ideas.',
        'The immutable v4 ledger and all original failures/source amendments are preserved in inputs.',
        'Frozen comparators are rerun on copied raw data; pilot/partial/missing evidence cannot establish PASS.',
        'Metric-specific errors retain their original units and meaning; A11 retains its L1 spectrum error.',
        'ED, thermal ED, DMRG and DQMC scope is explicit per row. A finite fixture does not certify all parameters.','',
        '| ID | Status | Scope or reason |','|---|---|---|']
    for row in rows:
        if row.get('input_role') in supplied and supplied[row['input_role']]:
            text.append(f"| {row['id']} | {row['status']} | {row['reason'] or row['executed_fixture_scope']} |")
    (output/'README.md').write_text('\n'.join(text)+'\n')
    print(json.dumps({'output_dir':str(output),'counts':counts,'audits':audits},indent=2))


if __name__=='__main__':
    main()
