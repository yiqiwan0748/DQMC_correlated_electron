#!/usr/bin/env python3
"""Independent arithmetic/provenance audit of the frozen L24 feasibility pilot.

No physics, submission or scientific catalogue completion. Stdlib only.
"""
import hashlib, json, math, re, shutil, tarfile, tomllib
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
RAW = BASE / 'results/bethe_l24_pilot_6458547'
OUT = BASE / 'reporting/reports/bethe_l24_audit_6458547'
CFG = BASE / 'bethe_ansatz_l24/pilot.toml'
checks = 0
def require(ok, label):
    global checks
    checks += 1
    if not ok:
        raise ValueError(label)
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def finite(x):
    if isinstance(x, dict): return all(finite(y) for y in x.values())
    if isinstance(x, list): return all(finite(y) for y in x)
    return not isinstance(x, float) or math.isfinite(x)

def main():
    require(not OUT.exists(), 'Use a new audit output')
    cfg = tomllib.loads(CFG.read_text())
    data = tomllib.loads((RAW/'target.toml').read_text())
    summary = json.loads((RAW/'summary.json').read_text())
    require(finite(data) and finite(summary), 'Finite receipts')
    require(cfg['L'] == [24] and cfg['scientific_completion_allowed'] is False, 'L24 pilot only')
    require(data['manifest_sha256'] == sha(CFG), 'Exact manifest')
    require(data['adapter_sha256'] == sha(BASE/'bethe_ansatz_l24/pilot_target.jl'), 'Actual adapter')
    require(data['status'] == 'PILOT_ONLY_COMPLETE', 'Completed pilot')
    for field, rel in [('source_sha256','analysis/compute_s2_itensor_dmrg.jl'),
                       ('project_sha256','Project.toml'), ('manifest_sha256','Manifest.toml')]:
        require(sha(BASE.parent/'validated_package/HubbardRenyiDQMC'/rel) == cfg[field], 'User source '+rel)
    freeze = json.loads((BASE/'source_releases/bethe-l24-pilot-v2-freeze.json').read_text())
    archive = BASE/'source_releases/bethe-l24-pilot-v2-source.tar.gz'
    require(sha(archive) == freeze['archive_sha256'], 'Source archive')
    seen = {}
    with tarfile.open(archive) as tf:
        for member in tf:
            require(member.isfile() and member.name not in seen, 'Regular unique frozen source')
            value = hashlib.sha256(tf.extractfile(member).read()).hexdigest()
            require(freeze['files'].get(member.name) == value, 'Frozen source '+member.name)
            seen[member.name] = value
    require(seen == freeze['files'], 'Full frozen source closure')
    transfers = {}
    remote = Path('/tmp/l24_remote_hashes_6458547.txt').read_text()
    for line in remote.splitlines():
        found = re.fullmatch(r'([a-f0-9]{64})  (.+)', line)
        if not found: continue
        digest, name = found.groups()
        path = BASE/name if name.startswith('logs/') else RAW/name
        require(sha(path) == digest, 'Remote transfer '+name)
        transfers[name] = digest
    require('6458547|COMPLETED|0:0|00:16:31|1|3500M' in remote, 'Scheduler identity/resources')
    rows = {(r['sector'],r['arm']):r for r in data['rows']}
    require(len(rows) == len(data['rows']) == 6, 'Six unique solves')
    require(set(rows) == {(s,a) for s in cfg['sectors'] for a in cfg['arms']}, 'Complete sectors/arms')
    g = cfg['gates']
    for (sector, arm), row in rows.items():
        require(tomllib.loads((RAW/f'L24_{sector}_{arm}.toml').read_text()) == row, 'Individual receipt identity')
        ac = cfg['arms'][arm]
        nu, nd = 12+int(sector=='plus'), 12-int(sector=='minus')
        require(row['seed'] == ac['seed_base']+24000+10*(cfg['sectors'].index(sector)+1), 'Predeclared seed')
        require(abs(row['norm']-1) <= g['norm_abs'], 'Norm')
        require(abs(row['Nup']-nu) <= g['sector_abs'] and abs(row['Ndown']-nd) <= g['sector_abs'], 'Particle/spin sector')
        require(abs(row['energy']-row['reported_energy']) <= g['reported_contracted_energy_abs'], 'Contracted energy')
        require(g['variance_min'] <= row['variance'] and row['variance']/24 <= ac['variance_per_site_max'], 'Variance')
        require(abs(row['variance_per_site']-row['variance']/24) < 1e-12, 'Variance units')
        require(abs(row['energy_per_site']-row['energy']/24) < 1e-12, 'Energy units')
        require(len(row['sweep_truncerrors']) == ac['nsweeps'], 'All sweeps retained')
        require(row['final_truncerr'] == row['sweep_truncerrors'][-1] and row['final_truncerr'] <= ac['final_truncerr_max'], 'Truncation')
        for key, total in [('site_Nup',row['Nup']),('site_Ndown',row['Ndown']),('site_double_occupancy',24*row['double_occupancy_per_site'])]:
            require(len(row[key]) == 24 and abs(sum(row[key])-total) <= 1e-9, 'Site normalization '+key)
    values = {}
    for arm in cfg['arms']:
        h, minus, plus = (rows[s,arm] for s in ('half','minus','plus'))
        require(abs(plus['energy']-minus['energy']-cfg['U']) <= g['particle_hole_energy_abs'], 'Particle-hole identity')
        values[arm] = {'energy_per_site':h['energy']/24, 'charge_gap':plus['energy']+minus['energy']-2*h['energy'], 'double_occupancy':sum(h['site_double_occupancy'])/24}
    errors = {k:abs(values['coarse'][k]-values['refined'][k]) for k in values['coarse']}
    for key in errors:
        require(errors[key] <= g['coarse_refined_'+key+'_abs'], 'Refinement stability '+key)
    require(summary == json.loads(Path('/tmp/l24_frozen_summary_6458547.json').read_text()), 'Exact frozen summary replay')
    require(summary['pilot_quality'] == 'PILOT_QUALITY_PASS' and summary['distinct_ideas_passed'] == 0, 'Pilot quality; zero completions')
    resources = (RAW/'target.resources.txt').read_text()
    rss = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', resources).group(1))
    eligibility = {'job_elapsed_seconds':991, 'peak_rss_kib':rss, 'wall_at_most_30_minutes':991<=1800, 'peak_at_most_2300_MiB':rss<=2300*1024}
    require(eligibility['wall_at_most_30_minutes'] and eligibility['peak_at_most_2300_MiB'], 'Frozen engineering bridge criteria')
    OUT.mkdir()
    shutil.copyfile('/tmp/l24_frozen_summary_6458547.json', OUT/'frozen_summary_replay.json')
    (OUT/'remote_sha256_and_scheduler.txt').write_text(remote)
    report = {'status':'PASS','checks':checks,'scientific_completions':0,'catalogue_status':{'A18':'NOT_RUN','A19':'NOT_RUN','A20':'NOT_RUN'},'refined':values['refined'],'coarse_refined_errors':errors,'next_L32_engineering_eligibility':eligibility,'remote_transfer_hashes':transfers,'auditor_sha256':sha(__file__)}
    (OUT/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'README.md').write_text('# L24 feasibility audit\n\nSix native user-MPO/ITensor DMRG solves passed the unchanged convergence gates. The exact frozen comparator replay and remote receipt hashes agree. Job wall time16m31s and process peakRSS1,484,548KiB satisfy the predeclared L32 engineering criteria. This permits a separately frozen larger pilot, not a thermodynamic scientific pass. A18/A19/A20 remain NOT_RUN.\n\nAll six sector/arm receipts and complete site/truncation arrays remain in the raw result directory. Coarse/refined charge gaps differ by about1.87e-7; the refined finite-size gap1.45911 still differs substantially from the infinite-chain oracle. Larger-size extrapolation and held-out validation remain required.\n')
    print(json.dumps({'status':'PASS','checks':checks,'errors':errors,'L32_eligible':eligibility,'scientific_completions':0}))

if __name__ == '__main__': main()
