"""Strict JSON metadata export; no numerical imports or solver work.

Only the pinned TeNPy initial Delta_E NaN receives a lossless diagnostic tag.
Pinned Config conversion is confined to the separate effective-parameters record.
Unsupported objects, ambiguous keys and cycles fail before JSON serialization.
"""
import json
import math
import os
from pathlib import Path
import struct
import sys

TAG_KIND = 'TENPY_INITIAL_UNDEFINED_DELTA_E_BINARY64'
SOURCE_SHA256 = 'b9ec50d49bd49d2e374878d37ce358b805cc48575f64aa5b8f3e9de7b9308e32'
CONFIG_SOURCE_SHA256 = '3f2b1e175f6974405b14f919a2e9a7a36edea2d1831f9747c8bb9ead6bf61549'
ALLOWED_PATH = ('sweep_statistics', 'Delta_E', 0)


def float64_bits(value):
    return struct.pack('>d', value).hex()


class ExportError(ValueError):
    pass


class NonFiniteExport(ExportError):
    def __init__(self, path, value):
        self.diagnostic = {'schema': 'P17_NONFINITE_EXPORT_DIAGNOSTIC_V1',
                           'status': 'UNRESOLVED', 'path': list(path),
                           'float64_bits_big_endian': float64_bits(value),
                           'classification': 'nan' if math.isnan(value) else 'infinity',
                           'scientific_credit': 0}
        super().__init__('Unregistered nonfinite export at ' + repr(tuple(path)))


class UnsupportedExport(ExportError):
    def __init__(self, path, value, reason='UNREGISTERED_TYPE'):
        typ = type(value)
        self.diagnostic = {'schema': 'P17_UNSUPPORTED_EXPORT_DIAGNOSTIC_V1',
                           'status': 'UNRESOLVED', 'path': list(path),
                           'type_module': typ.__module__,
                           'type_name': typ.__qualname__,
                           'reason': reason, 'scientific_credit': 0}
        super().__init__('Unsupported export at ' + repr(tuple(path)) + ': ' + reason)


def _keys(value, path):
    """Keep ordinary JSON layout, while rejecting coercion collisions."""
    rendered = set()
    for key in value:
        if type(key) not in (str, int):
            raise UnsupportedExport(path + ('<dictionary-key>',), key, 'UNREGISTERED_DICTIONARY_KEY_TYPE')
        name = str(key)
        if name in rendered:
            raise UnsupportedExport(path, value, 'AMBIGUOUS_JSON_DICTIONARY_KEYS')
        rendered.add(name)


def _walk(value, path, allow_initial, tagged, active=None):
    if active is None:
        active = set()
    if type(value) is float and not math.isfinite(value):
        if allow_initial and path == ALLOWED_PATH and math.isnan(value):
            tag = {'kind': TAG_KIND,
                   'float64_bits_big_endian': float64_bits(value),
                   'source_sha256': SOURCE_SHA256,
                   'reason': 'Pinned TeNPy first iteration sets E_old=np.nan; Delta_E=(E-E_old)/N_sweeps_check.'}
            tagged.append({'path': list(path), **tag})
            return tag
        raise NonFiniteExport(path, value)
    if value is None or type(value) in (str, bool, int, float):
        return value
    if type(value) not in (dict, list, tuple):
        raise UnsupportedExport(path, value)
    ident = id(value)
    if ident in active:
        raise UnsupportedExport(path, value, 'CYCLIC_CONTAINER')
    active.add(ident)
    try:
        if type(value) is dict:
            _keys(value, path)
            return {key: _walk(item, path + (key,), allow_initial, tagged, active)
                    for key, item in value.items()}
        items = [_walk(item, path + (i,), allow_initial, tagged, active)
                 for i, item in enumerate(value)]
        return tuple(items) if type(value) is tuple else items
    finally:
        active.remove(ident)


def prepare_sample(sample):
    """Copy a sample; preserve physics and tag only the known initial diagnostic."""
    tagged = []
    result = _walk(sample, (), True, tagged)
    if tagged:
        result['diagnostic_export_representation'] = {
            'schema': 'P17_TENPY_DIAGNOSTIC_NAN_TAG_V1',
            'tagged_entries': tagged,
            'physical_values_changed': False}
    return result


def _binary64_scalar(value, scalar_type, dtype, path):
    """Normalize ONLY an explicitly registered exact scalar class in metadata.

    Caller supplies the actually qualified numpy.float64 class. Physical exports
    never call this function. Original native bytes, dtype and path are retained.
    """
    if scalar_type is None or type(value) is not scalar_type:
        return None
    native_dtype=('<' if sys.byteorder=='little' else '>')+'f8'
    if dtype!=native_dtype:
        raise UnsupportedExport(path,value,'REGISTERED_BINARY64_DTYPE_NOT_NATIVE_F8')
    raw=scalar_type.tobytes(value)
    if type(raw) is not bytes or len(raw)!=8:
        raise UnsupportedExport(path,value,'REGISTERED_BINARY64_NATIVE_BYTES_INVALID')
    converted=float(value)
    if type(converted) is not float:
        raise UnsupportedExport(path,value,'REGISTERED_BINARY64_NOT_PYTHON_FLOAT')
    if not math.isfinite(converted):
        error=NonFiniteExport(path,converted)
        error.diagnostic.update(original_native_bytes_hex=raw.hex(),original_dtype=dtype)
        raise error
    if struct.pack('=d',converted)!=raw:
        raise UnsupportedExport(path,value,'REGISTERED_BINARY64_NATIVE_BITS_CHANGED')
    reloaded=json.loads(json.dumps(converted,allow_nan=False))
    if type(reloaded) is not float or struct.pack('=d',reloaded)!=raw:
        raise UnsupportedExport(path,value,'REGISTERED_BINARY64_JSON_BITS_CHANGED')
    return converted,{'path':list(path),'original_type_module':scalar_type.__module__,
                      'original_type_name':scalar_type.__qualname__,
                      'original_dtype':dtype,'original_native_bytes_hex':raw.hex(),
                      'converted_python_binary64_big_endian_hex':float64_bits(converted),
                      'native_bit_roundtrip_passed':True,'JSON_bit_roundtrip_passed':True}


def _config_metadata(value, config_type, path, configs, active, scalar_type, dtype):
    if type(value) is config_type:
        if id(value) in active:
            raise UnsupportedExport(path, value, 'CYCLIC_CONFIG')
        active.add(id(value))
        try:
            if type(value.name) is not str or type(value.unused) is not set or any(type(k) is not str for k in value.unused):
                raise UnsupportedExport(path, value, 'UNREGISTERED_CONFIG_METADATA')
            configs.append({'path': list(path), 'name': value.name,
                            'unused_keys': sorted(value.unused),
                            'source_sha256': CONFIG_SOURCE_SHA256})
            _config_metadata(value.options, config_type, path, configs, active, scalar_type, dtype)
        finally:
            active.remove(id(value))
    elif type(value) in (dict, list, tuple):
        if id(value) in active:
            raise UnsupportedExport(path, value, 'CYCLIC_CONTAINER')
        active.add(id(value))
        try:
            if type(value) is dict:
                _keys(value, path)
                for key, item in value.items():
                    _config_metadata(item, config_type, path + (key,), configs, active, scalar_type, dtype)
            else:
                for i, item in enumerate(value):
                    _config_metadata(item, config_type, path + (i,), configs, active, scalar_type, dtype)
        finally:
            active.remove(id(value))
    else:
        # Reject unknown objects and Config subclasses BEFORE pinned as_dict
        # can dispatch recursively through an unapproved subclass override.
        if _binary64_scalar(value,scalar_type,dtype,path) is None:
            _walk(value, path, False, [])


def _parameter_values(value, config_type, path, key_types, tuple_paths, scalar_type, dtype, scalar_records):
    if type(value) is config_type:
        # The actual caller binds Config.as_dict dispatch to the frozen source.
        value = config_type.as_dict(value)
    if type(value) is dict:
        _keys(value, path)
        if any(type(k) is int for k in value):
            key_types.append({'path': list(path), 'keys': [
                {'type': 'int' if type(k) is int else 'str', 'value': k} for k in value]})
        return {k: _parameter_values(v, config_type, path + (k,), key_types, tuple_paths, scalar_type, dtype, scalar_records)
                for k, v in value.items()}
    if type(value) in (list, tuple):
        if type(value) is tuple:
            tuple_paths.append(list(path))
        vals = [_parameter_values(v, config_type, path + (i,), key_types, tuple_paths, scalar_type, dtype, scalar_records)
                for i, v in enumerate(value)]
        return tuple(vals) if type(value) is tuple else vals
    normalized=_binary64_scalar(value,scalar_type,dtype,path)
    if normalized is not None:
        converted,record=normalized
        scalar_records.append(record)
        return converted
    return _walk(value, path, False, [])


def effective_parameters(value, config_type, *, binary64_type=None, binary64_dtype=None):
    """Return all post-call option values, Config metadata and lost JSON type metadata.

    No option is omitted. Input key types and tuple locations make JSON's ordinary
    key/sequence representation explicit. Unknown objects remain errors.
    """
    configs = []
    _config_metadata(value, config_type, (), configs, set(), binary64_type, binary64_dtype)
    key_types, tuple_paths = [], []
    scalar_records=[]
    vals = _parameter_values(value, config_type, (), key_types, tuple_paths, binary64_type, binary64_dtype, scalar_records)
    return {'schema': 'P17_TENPY_EFFECTIVE_PARAMETERS_V1', 'values': vals,
            'root_kind': 'actual_mutated_call_dictionary',
            'config_metadata': configs, 'mapping_key_types': key_types,
            'tuple_paths': tuple_paths, 'Config_source_sha256': CONFIG_SOURCE_SHA256,
            'effective_binary64_metadata_normalizations':scalar_records}


def verify_requested_parameters(requested, effective, path=()):
    """Every explicitly supplied call-boundary value must survive unchanged."""
    if type(requested) is dict:
        if type(effective) is not dict:
            raise UnsupportedExport(path, effective, 'REQUESTED_PARAMETER_STRUCTURE_CHANGED')
        for key, value in requested.items():
            if key not in effective:
                raise UnsupportedExport(path + (key,), effective, 'REQUESTED_PARAMETER_MISSING')
            verify_requested_parameters(value, effective[key], path + (key,))
    elif type(requested) in (list, tuple):
        if type(effective) is not type(requested) or len(requested) != len(effective):
            raise UnsupportedExport(path, effective, 'REQUESTED_PARAMETER_SEQUENCE_CHANGED')
        for i, value in enumerate(requested):
            verify_requested_parameters(value, effective[i], path + (i,))
    elif type(requested) is not type(effective) or requested != effective:
        raise UnsupportedExport(path, effective, 'REQUESTED_PARAMETER_VALUE_CHANGED')


def strict_json_text(value):
    """No nonfinite or unsupported object is admitted to persisted JSON."""
    checked = _walk(value, (), False, [])
    return json.dumps(checked, indent=2, allow_nan=False) + '\n'


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
