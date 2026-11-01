"""Strict JSON export with one lossless, source-established diagnostic NaN tag.

This module uses only the Python standard library and performs no numerical work.
It never admits nonfinite physical values. The only tagged location is the first
Delta_E entry in a TeNPy sweep-statistics table: the pinned 1.1.1 source forms it
from E_old=np.nan on the initial iteration. Actual NaN bits are retained.
"""
import json
import math
import os
from pathlib import Path
import struct

TAG_KIND = 'TENPY_INITIAL_UNDEFINED_DELTA_E_BINARY64'
SOURCE_SHA256 = 'b9ec50d49bd49d2e374878d37ce358b805cc48575f64aa5b8f3e9de7b9308e32'
ALLOWED_PATH = ('sweep_statistics', 'Delta_E', 0)


def float64_bits(value):
    return struct.pack('>d', value).hex()


class NonFiniteExport(ValueError):
    def __init__(self, path, value):
        self.diagnostic = {'schema': 'P17_NONFINITE_EXPORT_DIAGNOSTIC_V1',
                           'status': 'UNRESOLVED', 'path': list(path),
                           'float64_bits_big_endian': float64_bits(value),
                           'classification': 'nan' if math.isnan(value) else 'infinity',
                           'scientific_credit': 0}
        super().__init__('Unregistered nonfinite export at ' + repr(tuple(path)))


def _walk(value, path, allow_initial, tagged):
    if isinstance(value, float) and not math.isfinite(value):
        if allow_initial and path == ALLOWED_PATH and math.isnan(value):
            tag = {'kind': TAG_KIND,
                   'float64_bits_big_endian': float64_bits(value),
                   'source_sha256': SOURCE_SHA256,
                   'reason': 'Pinned TeNPy first iteration sets E_old=np.nan; Delta_E=(E-E_old)/N_sweeps_check.'}
            tagged.append({'path': list(path), **tag})
            return tag
        raise NonFiniteExport(path, value)
    if isinstance(value, dict):
        return {key: _walk(item, path + (key,), allow_initial, tagged)
                for key, item in value.items()}
    if isinstance(value, list):
        return [_walk(item, path + (i,), allow_initial, tagged)
                for i, item in enumerate(value)]
    if isinstance(value, tuple):
        return tuple(_walk(item, path + (i,), allow_initial, tagged)
                     for i, item in enumerate(value))
    return value


def prepare_sample(sample):
    """Return a copied sample with only the allowed metadata NaN tagged.

    Unknown/later NaNs and every infinity fail with an exact private path and
    original binary64 bits. Finite values are copied without conversion.
    """
    tagged = []
    result = _walk(sample, (), True, tagged)
    if tagged:
        result['diagnostic_export_representation'] = {
            'schema': 'P17_TENPY_DIAGNOSTIC_NAN_TAG_V1',
            'tagged_entries': tagged,
            'physical_values_changed': False}
    return result


def strict_json_text(value):
    """No nonfinite float is admitted anywhere in the persisted JSON value."""
    _walk(value, (), False, [])
    return json.dumps(value, indent=2, allow_nan=False) + '\n'


def write_new_json(path, value):
    text = strict_json_text(value)
    with Path(path).open('x', encoding='utf-8') as handle:
        handle.write(text)


def replace_checkpoint(path, value):
    """Atomically update only this new attempt's dedicated partial checkpoint."""
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    text = strict_json_text(value)
    with temporary.open('x', encoding='utf-8') as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
