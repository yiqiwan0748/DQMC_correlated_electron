#!/usr/bin/env python3
"""Meaningful standard-library-only strict-export boundary fixtures; no solver imports."""
import copy
import json
from pathlib import Path
import struct
import tempfile

from strict_export import (NonFiniteExport, TAG_KIND, float64_bits, prepare_sample,
                           strict_json_text, replace_checkpoint, write_new_json)


def fixtures():
    checks = []
    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)
    finite = {'fresh_energy': {'real': [-2.0], 'imag': [0.0]},
              'state': {'real': [0.0, 1.0]},
              'sweep_statistics': {'Delta_E': [0.0, -0.2], 'E': [-1.8, -2.0]}}
    control = copy.deepcopy(finite)
    check('finite sample remains identical', prepare_sample(finite) == control)
    check('input sample is not modified', finite == control)
    for bits in ['7ff8000000000000', 'fff800000000002a', '7ff8000000001234']:
        sample = copy.deepcopy(finite)
        value = struct.unpack('>d', bytes.fromhex(bits))[0]
        sample['sweep_statistics']['Delta_E'][0] = value
        actual = prepare_sample(sample)
        tag = actual['sweep_statistics']['Delta_E'][0]
        check('allowed initial NaN tag kind ' + bits, tag['kind'] == TAG_KIND)
        check('original NaN payload/sign bits retained ' + bits, tag['float64_bits_big_endian'] == bits)
        check('physical values exact ' + bits, actual['fresh_energy'] == finite['fresh_energy'] and actual['state'] == finite['state'])
        check('tag is strict RFC JSON ' + bits, json.loads(strict_json_text(actual)) == actual)
        check('original NaN was not overwritten ' + bits, float64_bits(sample['sweep_statistics']['Delta_E'][0]) == bits)
    failures = [
        ('later Delta_E NaN', ('sweep_statistics', 'Delta_E', 1), float('nan')),
        ('first Delta_E +infinity', ('sweep_statistics', 'Delta_E', 0), float('inf')),
        ('first Delta_E -infinity', ('sweep_statistics', 'Delta_E', 0), -float('inf')),
        ('physical energy NaN', ('fresh_energy', 'real', 0), float('nan')),
        ('physical state NaN', ('state', 'real', 1), float('nan')),
        ('other sweep-statistics NaN', ('sweep_statistics', 'E', 0), float('nan'))]
    for name, path, value in failures:
        sample = copy.deepcopy(finite); target = sample
        for part in path[:-1]: target = target[part]
        target[path[-1]] = value
        try:
            prepare_sample(sample)
        except NonFiniteExport as error:
            check(name + ' rejected', error.diagnostic['path'] == list(path))
            check(name + ' original bits retained', error.diagnostic['float64_bits_big_endian'] == float64_bits(value))
            check(name + ' diagnostic strict JSON', bool(json.loads(strict_json_text(error.diagnostic))))
        else:
            raise AssertionError(name + ' incorrectly admitted')
    try:
        strict_json_text({'physical': float('nan')})
    except NonFiniteExport:
        check('final report writer rejects unregistered NaN', True)
    else:
        raise AssertionError('writer admitted NaN')
    with tempfile.TemporaryDirectory() as folder:
        p = Path(folder) / 'partial.json'
        replace_checkpoint(p, {'status': 'RUNNING', 'samples': [finite]})
        before = p.read_bytes()
        try:
            replace_checkpoint(p, {'physical': float('nan')})
        except NonFiniteExport:
            check('failed next checkpoint preserves prior bytes', p.read_bytes() == before)
        else:
            raise AssertionError('nonfinite checkpoint admitted')
        replace_checkpoint(p, {'status': 'RUNNING', 'samples': [finite, finite]})
        check('next finite checkpoint atomically replaces previous', len(json.loads(p.read_text())['samples']) == 2)
        final = Path(folder) / 'final.json';write_new_json(final, finite)
        try:
            write_new_json(final, finite)
        except FileExistsError:
            check('final existing result never overwritten', True)
        else:
            raise AssertionError('existing result overwritten')
    return {'schema': 'P17_STRICT_EXPORT_METADATA_FIXTURES_V1', 'status': 'PASS',
            'checks': checks, 'check_count': len(checks), 'scientific_credit': 0,
            'solver_or_numerical_imports': False}

if __name__ == '__main__':
    print(json.dumps(fixtures(), indent=2, allow_nan=False))
