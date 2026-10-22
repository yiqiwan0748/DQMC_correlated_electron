#!/usr/bin/env python3
"""Verify this preserved research record without running a physics simulation."""
import csv
import hashlib
import json
from pathlib import Path
import tomllib

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'source_snapshot/HubbardRenyiDQMC'
ALM = ROOT / 'evidence/AltermagnetDQMCTrial'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rows(path):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))


def main():
    provenance = json.loads((ROOT / 'provenance.json').read_text())
    for entry in provenance['files']:
        file = ROOT / entry['record_path']
        require(file.is_file(), f'Missing preserved file: {file}')
        require(hashlib.sha256(file.read_bytes()).hexdigest() == entry['sha256'],
                f'SHA-256 mismatch: {file}')

    benchmark = rows(SOURCE / 'results/processed/BLOCKSIZE_U2_beta40_Ninc16_summary.csv')
    require([int(row['nA']) for row in benchmark] == [1, 2, 3, 4, 5], 'Subsystem coverage')
    zs = [(float(row['S2_DQMC']) - float(row['S2_ED'])) / float(row['S2_error'])
          for row in benchmark]
    require(max(map(abs, zs)) < 0.6, 'ED comparison no longer within 0.6 MC sigma')
    for row, z in zip(benchmark, zs):
        require(abs(z - float(row['z'])) < 0.00051, 'Printed z rounding mismatch')
        require(abs(float(row['DeltaS2']) - (float(row['S2_DQMC']) - float(row['S2_ED']))) < 2e-12,
                'Printed ED difference mismatch')

    dmrg = {int(row['L']): row for row in rows(SOURCE / 'results/processed/DMRG_S2_U2_periodic.csv')}
    large_l = {}
    for size in (20, 22):
        summary = {row['quantity']: row for row in rows(
            SOURCE / f'results/processed/FSS_1D_U2_beta4L_L{size}_Ninc16.csv')}
        value = float(summary['S2_DQMC']['value'])
        error = float(summary['S2_DQMC']['error'])
        ref = float(dmrg[size]['S2_MPS'])
        require(error > 0, 'Nonpositive MC error')
        large_l[str(size)] = {'S2_DQMC': value, 'nominal_MC_error': error,
                              'S2_DMRG': ref, 'z_using_DQMC_error_only': (value-ref)/error}

    fit_input = rows(SOURCE / 'results/processed/oneD_hubbard_fss_fit_input_L4_to_L22_20260820.csv')
    fits = rows(SOURCE / 'results/processed/oneD_hubbard_fss_fit_summary_L4_to_L22_20260820.csv')
    require(len(fit_input) == 22, 'FSS input row count')
    require(len(fits) == 289, 'FSS fit row count')
    principal = [row for row in fits if row['dataset_variant'] == 'dqmc' and row['role'] == 'principal']
    require(len(principal) == 5, 'Principal DQMC model count')
    chi2 = [float(row['statistic_value']) for row in principal]
    require(min(chi2) > 3, 'Previously poor full-window fit became adequate; review claim')

    attempt = ALM / 'results/literature_transfer/alm_dmrg_benchmark_20260911/attempt01'
    with (attempt / 'summary.toml').open('rb') as handle:
        dmrg_summary = tomllib.load(handle)
    audit_file = ALM / 'results/literature_transfer/alm_method_benchmark_20260911/dmrg_return_audit.toml'
    with audit_file.open('rb') as handle:
        audit = tomllib.load(handle)
    require(dmrg_summary['status'] == 'PASS' and audit['audit_status'] == 'PASS', 'Saved DMRG status')
    require(dmrg_summary['completed_arms'] == 4 and audit['final_pass_count'] == 4, 'Saved DMRG arms')
    for key in ('max_final_energy_error', 'max_final_spin_matrix_error', 'max_final_S2_error'):
        require(dmrg_summary[key] == audit['final_extrema'][key], f'Saved independent extrema: {key}')
    terminal = json.loads((ALM / 'results/isaac_20261014/moreo_configurations/solve003_terminal_validation_20261015.json').read_text())
    require(terminal['scientific_state_accepted'] is True, 'Saved ED acceptance')
    require(terminal['global_ground_state_certified'] is False and terminal['ALM_order_proved'] is False,
            'ED limitations changed; review claim')
    report = {
        'status': 'PASS_FILE_INTEGRITY_AND_RECORDED_SCALAR_ARITHMETIC',
        'preserved_files_verified': len(provenance['files']),
        'simulations_run_by_this_verifier': 0,
        'subsystem_z_recomputed': zs,
        'max_abs_subsystem_z': max(map(abs, zs)),
        'large_L_comparisons': large_l,
        'FSS_input_rows': len(fit_input), 'FSS_fit_rows': len(fits),
        'principal_DQMC_chi2_per_dof': chi2,
        'saved_DMRG_final_arms_passed': audit['final_pass_count'],
        'saved_DMRG_max_energy_error': dmrg_summary['max_final_energy_error'],
        'limits': ['No historical timing repeated', 'No large ED or MPS vector revalidated',
                   'No full external large-L raw production corpus recovered',
                   'No personal contribution or appointment dates inferred']}
    (ROOT / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
