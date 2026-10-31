#!/usr/bin/env python3
"""Independent arithmetic/provenance audit of the frozen L24 feasibility pilot.

No physics, submission or scientific catalogue completion. Stdlib only.
"""
import argparse, hashlib, json, math, re, shutil, subprocess, sys, tarfile, tempfile, tomllib
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
RAW = BASE / 'results/bethe_l24_pilot_6458547'
OUT = BASE / 'reporting/reports/bethe_l24_audit_6458547_portable_v2'
CFG = BASE / 'bethe_ansatz_l24/pilot.toml'
FIRST_AUDITOR_SHA = '5c6321310db6d21055aa4a421b2987ee58ddee870046418b32e5367b06a95dc1'
ENGINEERING_ANCHOR_SHA = '0e88dd7b1f051db016fbca825c2e2191872f93bda9184980c46b261d359dcea6'
ENGINEERING_ANCHOR_COMMIT = '66950c791f7f939ca176e5427be62cd8dfca8348'
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
    global RAW, OUT
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', type=Path, default=RAW)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args()
    RAW, OUT = args.raw.resolve(), args.output.resolve()
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
    require(sha(BASE/'reporting/auditor_versions/audit_bethe_l24_v1_6458547.py') == FIRST_AUDITOR_SHA, 'First auditor preserved')
    remote_input = RAW/'remote_sha256_and_scheduler_6458547.txt'
    replay_input = RAW/'frozen_summary_replay_6458547.json'
    require(sha(remote_input) == '8f13f8927d49e918550360ded36db0840cadeabe2a283994662d210a5edac9e6', 'Preserved original remote input')
    require(sha(replay_input) == '82fcbce021bcac1be43afbe0d08348dd36a67d064819ae77fdf47cca83c31b58', 'Preserved original frozen replay')
    anchor_path = RAW/'preexecution_L32_engineering_anchor_66950c7.md'
    anchor = json.loads((RAW/'preexecution_L32_engineering_anchor_66950c7.json').read_text())
    anchor_bytes = anchor_path.read_bytes()
    require(sha(anchor_path) == ENGINEERING_ANCHOR_SHA == anchor['sha256'], 'Exact published engineering anchor')
    require(anchor['commit'] == ENGINEERING_ANCHOR_COMMIT and anchor['commit_created_at'] == '2026-10-22T19:12:58Z', 'Published anchor identity/date')
    require(hashlib.sha1(b'blob '+str(len(anchor_bytes)).encode()+b'\x00'+anchor_bytes).hexdigest() == anchor['git_blob_sha'] == '6f09e3c24b16752df0bf6577ae9d9c631dc20af7', 'Published Git blob identity')
    require('at most 30 minutes and peak RSS at most 2300 MiB' in anchor_bytes.decode(), 'Predeclared engineering thresholds')
    transfers = {}
    scheduler_rows = []
    remote = remote_input.read_text()
    for line in remote.splitlines():
        found = re.fullmatch(r'([a-f0-9]{64})  (.+)', line)
        if not found:
            if line.startswith('6458547|'): scheduler_rows.append(line.split('|'))
            continue
        digest, name = found.groups()
        path = BASE/name if name.startswith('logs/') else RAW/name
        require(sha(path) == digest, 'Remote transfer '+name)
        transfers[name] = digest
    expected_transfers = {f'L24_{sector}_{arm}.toml' for sector in cfg['sectors'] for arm in cfg['arms']}
    expected_transfers |= {'target.toml','oracle.json','summary.json','oracle.resources.txt','target.resources.txt','logs/bethe_l24_pilot_6458547.log'}
    require(set(transfers) == expected_transfers, 'All original12 remote receipt/log hashes retained')
    require(len(scheduler_rows) == 1, 'Unique scheduler receipt')
    job_id,state,exit_code,elapsed,cpus,memory = scheduler_rows[0]
    require((job_id,state,exit_code,cpus,memory) == ('6458547','COMPLETED','0:0','1','3500M'), 'Scheduler identity/resources')
    elapsed_parts = [int(x) for x in elapsed.split(':')]
    require(len(elapsed_parts) == 3, 'Scheduler elapsed format')
    elapsed_seconds = 3600*elapsed_parts[0]+60*elapsed_parts[1]+elapsed_parts[2]
    rows = {(r['sector'],r['arm']):r for r in data['rows']}
    require(len(rows) == len(data['rows']) == 6, 'Six unique solves')
    require(set(rows) == {(s,a) for s in cfg['sectors'] for a in cfg['arms']}, 'Complete sectors/arms')
    g = cfg['gates']
    for (sector, arm), row in rows.items():
        require(tomllib.loads((RAW/f'L24_{sector}_{arm}.toml').read_text()) == row, 'Individual receipt identity')
        ac = cfg['arms'][arm]
        nu, nd = 12+int(sector=='plus'), 12-int(sector=='minus')
        require(row['L'] == 24 and row['Nup_expected'] == nu and row['Ndown_expected'] == nd, 'Expected L24 sector')
        require(row['seed'] == ac['seed_base']+24000+10*(cfg['sectors'].index(sector)+1), 'Predeclared seed')
        require(abs(row['norm']-1) <= g['norm_abs'], 'Norm')
        require(abs(row['Nup']-nu) <= g['sector_abs'] and abs(row['Ndown']-nd) <= g['sector_abs'], 'Particle/spin sector')
        require(abs(row['energy']-row['reported_energy']) <= g['reported_contracted_energy_abs'], 'Contracted energy')
        require(g['variance_min'] <= row['variance'] and row['variance']/24 <= ac['variance_per_site_max'], 'Variance')
        require(abs(row['variance_per_site']-row['variance']/24) < 1e-12, 'Variance units')
        require(abs(row['energy_per_site']-row['energy']/24) < 1e-12, 'Energy units')
        require(len(row['sweep_truncerrors']) == ac['nsweeps'], 'All sweeps retained')
        require(row['final_truncerr'] == row['sweep_truncerrors'][-1] and 0 <= row['final_truncerr'] <= ac['final_truncerr_max'], 'Truncation')
        require(0 < row['max_bond'] <= ac['chi_max'], 'Retained bond cap')
        if sector == 'half':
            require(g['double_occupancy_min'] <= row['double_occupancy_per_site'] <= g['double_occupancy_max'], 'Physical half-filled double occupancy')
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
    require(summary == json.loads(replay_input.read_text()), 'Preserved exact frozen summary replay')
    # Re-execute only the frozen stdlib summarizer, never a physics driver.
    with tempfile.TemporaryDirectory(prefix='bethe_l24_postprocess_') as tempdir:
        temp_root = Path(tempdir)
        with tarfile.open(archive) as tf:
            for member in tf:
                rel = Path(member.name)
                require(not rel.is_absolute() and '..' not in rel.parts, 'Safe frozen source member')
                dest = temp_root/rel
                dest.parent.mkdir(parents=True,exist_ok=True)
                dest.write_bytes(tf.extractfile(member).read())
        frozen_root = temp_root/'bethe_ansatz_l24'
        replay_output = temp_root/'fresh_summary_replay.json'
        subprocess.run([sys.executable,str(frozen_root/'summarize_pilot.py'),
            '--manifest',str(frozen_root/'pilot.toml'),'--oracle',str(RAW/'oracle.json'),
            '--target',str(RAW/'target.toml'),'--output',str(replay_output)],
            check=True,capture_output=True,text=True)
        fresh_replay = replay_output.read_bytes()
    require(fresh_replay == replay_input.read_bytes() == (RAW/'summary.json').read_bytes(), 'Fresh frozen replay agrees byte-for-byte')
    require(summary['pilot_quality'] == 'PILOT_QUALITY_PASS' and summary['distinct_ideas_passed'] == 0, 'Pilot quality; zero completions')
    resources = (RAW/'target.resources.txt').read_text()
    rss = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', resources).group(1))
    require('Exit status: 0' in resources, 'Successful target resource receipt')
    eligibility = {'job_elapsed_seconds':elapsed_seconds, 'peak_rss_kib':rss, 'peak_rss_MiB':rss/1024,
                   'wall_at_most_30_minutes':elapsed_seconds<=1800, 'peak_at_most_2300_MiB':rss<=2300*1024,
                   'threshold_scope':'Predeclared engineering continuation criteria; separate from unchanged L24 scientific-quality gates',
                   'preexecution_anchor':anchor}
    require(eligibility['wall_at_most_30_minutes'] and eligibility['peak_at_most_2300_MiB'], 'Frozen engineering bridge criteria')
    OUT.mkdir()
    (OUT/'frozen_summary_replay.json').write_bytes(fresh_replay)
    shutil.copyfile(__file__,OUT/'auditor_source.py')
    (OUT/'remote_sha256_and_scheduler.txt').write_text(remote)
    sidecars = [remote_input,replay_input,anchor_path,RAW/'preexecution_L32_engineering_anchor_66950c7.json',RAW/'raw_sidecars_provenance_6458547.json']
    report = {'version':'portable_v2','status':'PASS','checks':checks,'scientific_completions':0,'catalogue_status':{'A18':'NOT_RUN','A19':'NOT_RUN','A20':'NOT_RUN'},'refined':values['refined'],'coarse_refined_errors':errors,'next_L32_engineering_eligibility':eligibility,'remote_transfer_hashes':transfers,'raw_sidecar_hashes':{p.name:sha(p) for p in sidecars},'auditor_sha256':sha(__file__),'first_auditor_sha256':FIRST_AUDITOR_SHA,'frozen_summary_replay_sha256':hashlib.sha256(fresh_replay).hexdigest(),'physics_executed':False}
    (OUT/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'README.md').write_text('# Portable L24 feasibility audit, revision2\n\nSix native user-MPO/ITensor DMRG solves pass the unchanged convergence gates. All12 original remote scientific/log hashes agree. A fresh replay of the exact summarizer extracted from the frozen source archive agrees byte-for-byte with both the original raw summary and preserved replay. This is arithmetic/postprocessing only; no numerical physics was rerun.\n\nThe original auditor and first report are preserved. Required original hash/scheduler and replay inputs are now named sidecars in the raw result directory; no temporary input path is required. The two historical input sidecars retain exactly their original bytes.\n\nThe30-minute/2300-MiB L32 engineering criteria are bound to the exact previously published README at commit66950c791f7f939ca176e5427be62cd8dfca8348 (2026-10-22T19:12:58Z), preserved as a raw sidecar with its Gitblob and SHA256. Job elapsed time16m31s=991s and target peakRSS1,484,548KiB=1449.75390625MiB satisfy those criteria. They are resource continuation gates, separate from scientific convergence tolerances.\n\nAll six sector/arm receipts and complete site/truncation arrays remain in the raw directory. Coarse/refined charge gaps differ by about1.87e-7; the refined finite-size gap1.45911 still differs substantially from the infinite-chain oracle. This permits a separately frozen larger feasibility pilot; it establishes no thermodynamic result. A18/A19/A20 remain NOT_RUN and scientific completions remain0.\n')
    print(json.dumps({'status':'PASS','checks':checks,'errors':errors,'L32_eligible':eligibility,'scientific_completions':0}))

if __name__ == '__main__': main()
