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



def config_and_type_fixtures():
    """Metadata-only stand-in fixtures; actual Config runtime remains unqualified here."""
    from strict_export import (ExportError, UnsupportedExport, effective_parameters,
                               verify_requested_parameters)
    checks=[]
    def check(name,value):
        if not value: raise AssertionError(name)
        checks.append(name)
    class ConfigStandIn:
        # Deliberate metadata fixture implementing the inspected Config value contract;
        # this is not the real package or a runtime-qualification claim.
        def __init__(self,options,name,unused=()):
            self.options=options;self.name=name;self.unused=set(unused)
        def as_dict(self):
            result=dict(self.options)
            for key,value in result.items():
                if isinstance(value,ConfigStandIn):result[key]=value.as_dict()
            return result
    requested={'min_sweeps':15,'max_sweeps':15,
               'chi_list':{0:8,1:16,2:32,3:64},
               'trunc_params':{'chi_max':64,'svd_min':1e-16,'trunc_cut':1e-8}}
    actual=copy.deepcopy(requested)
    actual['trunc_params']=ConfigStandIn({**actual['trunc_params'],'chi_min':None,
                                          'degeneracy_tol':None},'trunc_params',('chi_min',))
    actual['added_default']=True
    exported=effective_parameters(actual,ConfigStandIn)
    verify_requested_parameters(requested,exported['values'])
    check('all explicit requested leaves retained',True)
    check('effective defaults retained without filtering',exported['values']['trunc_params']['chi_min'] is None and exported['values']['trunc_params']['degeneracy_tol'] is None and exported['values']['added_default'] is True)
    check('Config name and unused metadata retained',exported['config_metadata'][0]['name']=='trunc_params' and exported['config_metadata'][0]['unused_keys']==['chi_min'])
    check('integer chi-list key types retained separately',exported['mapping_key_types']==[{'path':['chi_list'],'keys':[{'type':'int','value':i} for i in range(4)]}])
    text=strict_json_text(exported);reloaded=json.loads(text)
    check('ordinary JSON key representation remains comparator compatible',reloaded['values']['chi_list']=={'0':8,'1':16,'2':32,'3':64})
    check('requested dict not mutated by effective conversion',requested['trunc_params']=={'chi_max':64,'svd_min':1e-16,'trunc_cut':1e-8})
    check('actual Config not replaced or altered',type(actual['trunc_params']) is ConfigStandIn and actual['trunc_params'].unused=={'chi_min'})
    nested=ConfigStandIn({'child':ConfigStandIn({'value':3},'child')},'outer')
    nested_export=effective_parameters({'settings':nested},ConfigStandIn)
    check('all nested Config names retained',[x['name'] for x in nested_export['config_metadata']]==['outer','child'])
    check('nested Config values retained',nested_export['values']=={'settings':{'child':{'value':3}}})
    sequences=effective_parameters({'nested':[ConfigStandIn({'x':1},'list-config')], 'pair':(1,2)},ConfigStandIn)
    check('Config inside sequences converted explicitly',sequences['values']['nested']==[{'x':1}])
    check('tuple paths preserve JSON-lost sequence type',sequences['tuple_paths']==[['pair']])
    def must_fail(name,action,path,reason):
        try:action()
        except UnsupportedExport as error:
            check(name+' exact path',error.diagnostic['path']==path)
            check(name+' reason',error.diagnostic['reason']==reason)
            check(name+' primitive diagnostic',json.loads(strict_json_text(error.diagnostic))==error.diagnostic)
        else:raise AssertionError(name+' admitted')
    bad=object()
    must_fail('unknown physical object',lambda:prepare_sample({'state':{'real':[bad]}}),['state','real',0],'UNREGISTERED_TYPE')
    must_fail('unknown effective option',lambda:effective_parameters({'settings':bad},ConfigStandIn),['settings'],'UNREGISTERED_TYPE')
    class ConfigSubclass(ConfigStandIn):
        override_called=False
        def as_dict(self):
            ConfigSubclass.override_called=True
            return {'silently_normalized':'unapproved subclass'}
    nested_override=ConfigStandIn({'nested':ConfigSubclass({},'bad-subclass')},'outer')
    must_fail('Config subclass rejected before recursive as_dict',
              lambda:effective_parameters({'settings':nested_override},ConfigStandIn),
              ['settings','nested'],'UNREGISTERED_TYPE')
    check('unapproved subclass as_dict override never called',not ConfigSubclass.override_called)
    must_fail('unconverted Config cannot leak through strict writer',lambda:strict_json_text({'settings':actual['trunc_params']}),['settings'],'UNREGISTERED_TYPE')
    cyclic=[];cyclic.append(cyclic)
    must_fail('cyclic list rejected',lambda:strict_json_text({'cycle':cyclic}),['cycle',0],'CYCLIC_CONTAINER')
    cycle_config=ConfigStandIn({},'cycle');cycle_config.options['again']=cycle_config
    must_fail('cyclic Config rejected',lambda:effective_parameters({'settings':cycle_config},ConfigStandIn),['settings','again'],'CYCLIC_CONFIG')
    must_fail('ambiguous integer/string keys rejected',lambda:strict_json_text({1:'one','1':'string-one'}),[],'AMBIGUOUS_JSON_DICTIONARY_KEYS')
    must_fail('unsupported key type rejected',lambda:strict_json_text({None:1}),['<dictionary-key>'],'UNREGISTERED_DICTIONARY_KEY_TYPE')
    changed=copy.deepcopy(exported['values']);changed['trunc_params']['chi_max']=32
    must_fail('changed requested leaf rejected',lambda:verify_requested_parameters(requested,changed),['trunc_params','chi_max'],'REQUESTED_PARAMETER_VALUE_CHANGED')
    missing=copy.deepcopy(exported['values']);del missing['max_sweeps']
    must_fail('missing requested leaf rejected',lambda:verify_requested_parameters(requested,missing),['max_sweeps'],'REQUESTED_PARAMETER_MISSING')
    with tempfile.TemporaryDirectory() as folder:
        partial=Path(folder)/'partial.json';report={'status':'RUNNING','samples':[{'good':True}]}
        replace_checkpoint(partial,report);before=partial.read_bytes()
        try:
            invalid=prepare_sample({'actual_dmrg_parameters':bad})
            candidate={**report,'samples':[*report['samples'],invalid]}
            strict_json_text(candidate);replace_checkpoint(partial,candidate);report=candidate
        except ExportError as error:
            final={**report,'status':'UNRESOLVED','error_diagnostic':error.diagnostic}
            write_new_json(Path(folder)/'final.json',final)
        else:raise AssertionError('unknown option admitted')
        check('invalid next sample preserves previous partial',partial.read_bytes()==before)
        final=json.loads((Path(folder)/'final.json').read_text())
        check('invalid next sample preserves primitive final report',final['status']=='UNRESOLVED' and final['samples']==[{'good':True}])
        check('bad new sample is not committed',report['samples']==[{'good':True}])
    return checks


def binary64_metadata_fixtures():
    """Exact bit-contract stand-ins; no actual NumPy import or runtime claim."""
    import sys
    from strict_export import (UnsupportedExport, effective_parameters,
                               verify_requested_parameters)
    checks=[]
    def check(name,condition):
        if not condition:raise AssertionError(name)
        checks.append(name)
    class ConfigStandIn:
        def __init__(self,options):
            self.options=options;self.name='lanczos';self.unused=set()
        def as_dict(self):return dict(self.options)
    class Binary64StandIn:
        calls=0
        def __init__(self,big_endian_hex):
            raw=bytes.fromhex(big_endian_hex)
            self.raw=raw[::-1] if sys.byteorder=='little' else raw
        def tobytes(self):
            Binary64StandIn.calls+=1
            return self.raw
        def __float__(self):
            Binary64StandIn.calls+=1
            return struct.unpack('=d',self.raw)[0]
    dtype=('<' if sys.byteorder=='little' else '>')+'f8'
    patterns=['0000000000000000','8000000000000000','0000000000000001',
              '0010000000000000','3ff0000000000000','7fefffffffffffff',
              'ffefffffffffffff']
    for bits in patterns:
        scalar=Binary64StandIn(bits)
        root={'lanczos_params':ConfigStandIn({'cutoff':scalar,'N_max':40})}
        actual=effective_parameters(root,ConfigStandIn,
                    binary64_type=Binary64StandIn,binary64_dtype=dtype)
        value=actual['values']['lanczos_params']['cutoff']
        row=actual['effective_binary64_metadata_normalizations'][0]
        check('registered metadata normalized exact binary64 '+bits,float64_bits(value)==bits)
        check('original native bytes retained '+bits,row['original_native_bytes_hex']==scalar.raw.hex())
        check('original dtype and path retained '+bits,row['original_dtype']==dtype and row['path']==['lanczos_params','cutoff'])
        check('Python bits and both roundtrips retained '+bits,row['converted_python_binary64_big_endian_hex']==bits and row['native_bit_roundtrip_passed'] and row['JSON_bit_roundtrip_passed'])
        reloaded=json.loads(strict_json_text(actual))
        check('persisted JSON preserves scalar bits '+bits,float64_bits(reloaded['values']['lanczos_params']['cutoff'])==bits)
        verify_requested_parameters({'lanczos_params':{'N_max':40}},actual['values'])
    def fail(name,action,error_type,path,reason=None):
        try:action()
        except error_type as error:
            check(name+' exact path',error.diagnostic['path']==path)
            if reason is not None:check(name+' reason',error.diagnostic['reason']==reason)
            check(name+' primitive diagnostic',json.loads(strict_json_text(error.diagnostic))==error.diagnostic)
        else:raise AssertionError(name+' admitted')
    for bits in ['7ff800000000002a','fff0000000000000','7ff0000000000000']:
        scalar=Binary64StandIn(bits)
        try:effective_parameters({'cutoff':scalar},ConfigStandIn,
                  binary64_type=Binary64StandIn,binary64_dtype=dtype)
        except NonFiniteExport as error:
            check('nonfinite registered metadata rejected '+bits,error.diagnostic['path']==['cutoff'])
            check('nonfinite original native bits retained '+bits,error.diagnostic['original_native_bytes_hex']==scalar.raw.hex())
        else:raise AssertionError('nonfinite registered metadata admitted')
    class Binary64Subclass(Binary64StandIn):
        override_calls=0
        def tobytes(self):
            Binary64Subclass.override_calls+=1
            return super().tobytes()
        def __float__(self):
            Binary64Subclass.override_calls+=1
            return super().__float__()
    subclass=Binary64Subclass('3ff0000000000000')
    fail('scalar subclass rejected before conversion',
         lambda:effective_parameters({'cutoff':subclass},ConfigStandIn,
                     binary64_type=Binary64StandIn,binary64_dtype=dtype),
         UnsupportedExport,['cutoff'],'UNREGISTERED_TYPE')
    check('scalar subclass methods never called',Binary64Subclass.override_calls==0)
    original=Binary64StandIn('3ff0000000000000');before=Binary64StandIn.calls
    fail('registered scalar still rejected in physical export',
         lambda:prepare_sample({'fresh_energy':{'real':[original]}}),
         UnsupportedExport,['fresh_energy','real',0],'UNREGISTERED_TYPE')
    check('physical rejection never invokes scalar conversion',Binary64StandIn.calls==before)
    fail('missing scalar registration remains closed',
         lambda:effective_parameters({'cutoff':original},ConfigStandIn),
         UnsupportedExport,['cutoff'],'UNREGISTERED_TYPE')
    fail('non-native dtype registration rejected',
         lambda:effective_parameters({'cutoff':original},ConfigStandIn,
                   binary64_type=Binary64StandIn,binary64_dtype=('>f8' if dtype=='<f8' else '<f8')),
         UnsupportedExport,['cutoff'],'REGISTERED_BINARY64_DTYPE_NOT_NATIVE_F8')
    class WrongBytes(Binary64StandIn):
        def tobytes(self):return b'bad'
    fail('wrong native byte count rejected',
         lambda:effective_parameters({'cutoff':WrongBytes('3ff0000000000000')},ConfigStandIn,
                  binary64_type=WrongBytes,binary64_dtype=dtype),
         UnsupportedExport,['cutoff'],'REGISTERED_BINARY64_NATIVE_BYTES_INVALID')
    class WrongFloat(Binary64StandIn):
        def __float__(self):return 2.0
    fail('nonlossless registered conversion rejected',
         lambda:effective_parameters({'cutoff':WrongFloat('3ff0000000000000')},ConfigStandIn,
                  binary64_type=WrongFloat,binary64_dtype=dtype),
         UnsupportedExport,['cutoff'],'REGISTERED_BINARY64_NATIVE_BITS_CHANGED')
    return checks


if __name__ == '__main__':
    original=fixtures()
    original['checks'].extend(config_and_type_fixtures())
    original['checks'].extend(binary64_metadata_fixtures())
    original['schema']='P17_STRICT_EXPORT_METADATA_FIXTURES_V3'
    original['check_count']=len(original['checks'])
    original['Config_stand_in_is_not_actual_package_qualification']=True
    original['binary64_stand_in_is_not_actual_NumPy_qualification']=True
    print(json.dumps(original,indent=2,allow_nan=False))
