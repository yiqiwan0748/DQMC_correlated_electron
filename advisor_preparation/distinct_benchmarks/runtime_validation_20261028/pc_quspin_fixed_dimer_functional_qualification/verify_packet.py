"""Verify public bytes and saved-data projection; no physical calculation."""
import ast
import hashlib
import json
import math
from pathlib import Path
import sys
from fractions import Fraction


ROOT = Path(__file__).resolve().parent


def need(value, message):
    if not value:
        raise ValueError(message)


def read(name):
    def real(text):
        value = float(text)
        need(math.isfinite(value), 'Nonfinite JSON scalar')
        return value
    def pairs(items):
        result = {}
        for key, value in items:
            need(key not in result, 'Duplicate JSON key')
            result[key] = value
        return result
    return json.loads((ROOT / name).read_bytes(), parse_float=real, object_pairs_hook=pairs,
        parse_constant=lambda text: (_ for _ in ()).throw(ValueError('Nonfinite JSON')))


def pin(name):
    raw = (ROOT / name).read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def rational(record):
    need(type(record) is dict and set(record) == {'numerator', 'denominator'}, 'Rational record required')
    n, d = record['numerator'], record['denominator']
    need(type(n) is type(d) is str and str(int(n)) == n and str(int(d)) == d and int(d) > 0,
         'Canonical integer rational required')
    value = Fraction(int(n), int(d))
    need((str(value.numerator), str(value.denominator)) == (n, d), 'Reduced rational required')
    return value


def vector(value, size):
    return type(value) is list and len(value) == size and all(type(x) is float and math.isfinite(x) for x in value)


def matrix(value, size):
    return type(value) is list and len(value) == size and all(vector(row, size) for row in value)


def main():
    manifest = read('MANIFEST.json')
    files = manifest['files']
    need(not any(path.is_symlink() for path in ROOT.rglob('*')), 'Packet symlink')
    need({path.relative_to(ROOT).as_posix() for path in ROOT.rglob('*') if path.is_file()} ==
         set(files) | {'MANIFEST.json'}, 'Packet file census differs')
    need(all(pin(name) == descriptor for name, descriptor in files.items()), 'Packet bytes differ')
    projection = read('PROJECTION_AUDIT.json')
    saved = read('saved_scientific_record.json')
    data = read('scientific_data.json')
    protocol = read('DIMER_PROTOCOL.json')
    need(all(data[key] == saved[key] for key in projection['unchanged_scientific_fields']), 'Saved projection differs')
    need(all(pin(name) == descriptor for name, descriptor in projection['scientific_file_pins'].items()),
         'Scientific source/checkpoint identity differs')
    need(data['model'] == protocol['model'] and data['scientific_credit'] == 0 and data['new_ideas'] == 0 and
         data['catalogue_count_changed'] is False and data['X08_status'] == 'NOT_RUN', 'Scope differs')
    need(vector(data['E'], 4) and matrix(data['H'], 4) and matrix(data['V_by_column'], 4) and
         len(data['actual_number_operators']) == 4 and all(matrix(item['matrix'], 4) for item in data['actual_number_operators']),
         'Incomplete or nonfinite saved arrays')
    need(len(data['analytic_energy_gates']) == 4 and all(item['index'] == i and item['pass'] is True and
         rational(item['limit']) == Fraction(1, 10**11) and
         0 <= rational(item['exact_distance_upper']) <= rational(item['limit'])
         for i, item in enumerate(data['analytic_energy_gates'])), 'Recorded energy bound fails')
    need(rational(data['Hermitian_limit']) == Fraction(1, 10**12) and
         0 <= rational(data['Hermitian_max_abs']) <= rational(data['Hermitian_limit']), 'Recorded Hermitian bound fails')
    need(data['eigenvector_global_gap_ground_state_qualification'] == 'NOT_CERTIFIED' and
         data['diagnostics']['qualification'] == 'POINT_DIAGNOSTICS_NOT_EIGENVECTOR_OR_GLOBAL_GAP_CERTIFICATES',
         'Diagnostic scope differs')
    ast.parse((ROOT / 'dimer_science.py').read_text())
    need(not {'numpy', 'scipy', 'quspin', 'quspin_extensions'}.intersection(sys.modules), 'Unexpected scientific import')
    print(json.dumps({'status': 'PASS_PUBLIC_BYTES_AND_SAVED_PROJECTION', 'files': len(files) + 1,
                      'physical_computations': 0, 'new_ideas': 0, 'full_X08': 'NOT_RUN'}))


if __name__ == '__main__':
    main()
