#!/usr/bin/env python3
"""Full integer identity comparisons of saved diagnostic metadata only."""
import argparse, json, tomllib
from pathlib import Path

def main():
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);x=a.parse_args();root=Path(__file__).resolve().parent
    j=tomllib.loads((x.output/'julia_observations_PRIVATE.toml').read_text());p=json.loads((x.output/'python_C_GNU_observations_PRIVATE.json').read_text())
    before={v['path']:v for v in j.get('before_Base_stat_rows',[])};after={v['path']:v for v in j.get('after_Base_stat_rows',[])}
    report=dict(schema='P17_FULLWIDTH_STAT_METADATA_COMPARISONS_V1',status='OBSERVATIONS_COMPLETE_NOT_RUNTIME_QUALIFICATION',scientific_completions=0,runtime_qualified=False,prospective_pin_approval=False,rows=[],current_mapped_segments=[],historical_v2_segments=[],errors=[])
    GNU_ok=all(p.get(k,{}).get('complete_exact_order_and_decimal_strings')is True for k in ['GNU_before','GNU_after'])
    paths=j['paths'];GNU={q:(p['GNU_before']['rows'][i],p['GNU_after']['rows'][i])for i,q in enumerate(paths)}if GNU_ok else{}
    for row in p['rows']:
        q=row['path'];b=before.get(q,{});e=after.get(q,{});out=dict(path=q,realpath=row.get('realpath'),observation=row,Base_before=b,Base_after=e);report['rows'].append(out)
        try:
            py=row['python_before'];py2=row['python_after'];c=row['C_before']['result'];c2=row['C_after']['result'];g,g2=GNU[q]
            full={k:int(py[k])for k in ['inode','device','major','minor','bytes']}
            if not all(isinstance(py[k],str)and py[k].isdecimal()and 0<=full[k]<2**64 for k in full):raise ValueError('Not full exact unsigned decimal strings')
            c_match=all(row[z]['exit_code']==0 for z in ['C_before','C_after'])and all(c[k+'_decimal']==c2[k+'_decimal']==py[k]for k in full)and c['fd_inode_decimal']==c['inode_decimal']and c['fd_device_decimal']==c['device_decimal']and c2['fd_inode_decimal']==c2['inode_decimal']and c2['fd_device_decimal']==c2['device_decimal']
            g_match=g==g2==[py['device'],py['inode'],py['bytes']]
            stable=py==py2 and row['python_before_after_equal']and b.get('sha256')==b.get('sha256_after')==e.get('sha256')==e.get('sha256_after')==row['sha256']
            stable=stable and b.get('Base_before_after_identity_equal')is True and e.get('Base_before_after_identity_equal')is True
            for field,tail in [('inode','inode'),('device','device'),('size','size')]:
                stable=stable and b.get('julia_'+field+'_value_string')==b.get('julia_after_'+tail+'_string')==e.get('julia_'+field+'_value_string')==e.get('julia_after_'+tail+'_string')
            stable=stable and all(b.get('julia_'+field+'_value_string')==e.get('julia_'+field+'_value_string')for field in ['mtime','ctime'])
            out.update(full_integer_C_Python_GNU_equal=c_match and g_match,all_hash_and_before_after_states_equal=stable,full_values=py,
                Base_full_inode_equal=int(b['julia_inode_value_string'])==full['inode'],Base_full_device_equal=int(b['julia_device_value_string'])==full['device'],
                Base_inode_is_full_mod_2pow32_DIAGNOSTIC_ONLY=int(b['julia_inode_value_string'])==full['inode']%(2**32),
                Base_device_is_full_mod_2pow32_DIAGNOSTIC_ONLY=int(b['julia_device_value_string'])==full['device']%(2**32))
        except Exception as err:out['comparison_error']=type(err).__name__+': '+str(err)
    bypath={v['path']:v for v in report['rows']}
    def segment_result(s,historical=False):
        q=s['path'];out=dict(segment=s,historical_metadata_only=historical)
        if not q.startswith('/'):return None
        v=bypath.get(q)
        if v is None:
            v=next((z for z in report['rows']if z.get('realpath')==q),None)
        if v and 'full_values'in v:
            f=v['full_values'];out['full_C_Python_GNU_identity_matches_map']=int(s['inode'])==int(f['inode'])and[int(z,16)for z in s['device'].split(':')]==[int(f['major']),int(f['minor'])]
            out['sha256']=v['observation'].get('sha256');out['all_stat_arms_equal']=v['full_integer_C_Python_GNU_equal'];out['before_after_states_equal']=v['all_hash_and_before_after_states_equal']
        else:out['unresolved']='Missing complete independent full-width stat values'
        return out
    for phase in ['before','after']:
        for s in j.get(phase+'_segments',[]):
            v=segment_result(s)
            if v:report['current_mapped_segments'].append(dict(phase=phase,**v))
    old=tomllib.loads((root/'metadata/v2_report_PRIVATE.toml').read_text())
    for s in old['mappings']:
        v=segment_result(s,True)
        if v:report['historical_v2_segments'].append(v)
    report['same_process_before_after_file_census_equal']=j.get('before_after_file_census_equal',False)
    report['all_current_full_identity_checks_agree']=bool(report['current_mapped_segments'])and all(v.get('full_C_Python_GNU_identity_matches_map')and v.get('all_stat_arms_equal')and v.get('before_after_states_equal')for v in report['current_mapped_segments'])
    historical=json.loads((root/'metadata/historical_paths.json').read_text())
    NAS=[v for v in report['rows']if v['path']in historical['unique_paths']and v['path'].startswith('/lustre/')]
    unknown=[v for v in report['rows']if 'Unregistered provider'in historical['classes'].get(v['path'],{}).get('historical_errors',[])]
    report['historical_NAS_paths_count']=len(NAS)
    report['all_current_stat_arms_for_historical_NAS_paths_equal_and_stable']=len(NAS)==24 and all(v.get('full_integer_C_Python_GNU_equal')and v.get('all_hash_and_before_after_states_equal')for v in NAS)
    report['NAS_full_inode_bit_lengths']={v['path']:int(v['full_values']['inode']).bit_length()for v in NAS if'full_values'in v}
    report['historical_unknown_provider_path_count']=len(unknown)
    report['historical_unknown_provider_current_observations']=[{'path':v['path'],'sha256':v['observation'].get('sha256'),'full_C_Python_GNU_equal':v.get('full_integer_C_Python_GNU_equal'),'prospective_pin_approval':False}for v in unknown]
    report['historical_unknown_paths_are_not_approved']=True
    report['full_value_disagreement_or_unstable_paths']=[v['path']for v in report['rows']if not(v.get('full_integer_C_Python_GNU_equal')and v.get('all_hash_and_before_after_states_equal'))]
    required=dict(observation_sources_complete=j.get('status')=='OBSERVATIONS_COMPLETE_NOT_RUNTIME_QUALIFICATION'and p.get('status')=='OBSERVATIONS_COMPLETE',
        exact_registered_path_census=bool(paths)and len(paths)==len(set(paths))==len(report['rows'])and set(paths)==set(bypath)==set(before)==set(after),
        every_historical_registered_path_sampled=set(historical['unique_paths']).issubset(set(paths)),
        both_GNU_batches_complete=GNU_ok,all_registered_path_full_arms_equal_and_stable=bool(report['rows'])and not report['full_value_disagreement_or_unstable_paths']and all('comparison_error'not in v for v in report['rows']),
        both_current_segment_phases_nonempty=all(any(v['phase']==phase for v in report['current_mapped_segments'])for phase in ['before','after']),
        current_full_identity_checks_agree=report['all_current_full_identity_checks_agree'],
        before_after_current_file_census_equal=report['same_process_before_after_file_census_equal'],
        C_helper_before_after_bytes_equal=isinstance(j.get('compiled_C_helper_sha256'),str)and len(j['compiled_C_helper_sha256'])==64 and all(v in'0123456789abcdef'for v in j['compiled_C_helper_sha256'])and j['compiled_C_helper_sha256']==p.get('helper_sha256')==p.get('helper_after_sha256'))
    report['required_metadata_closure_checks']=required
    report['required_metadata_closure']=all(required.values())
    report['status']='METADATA_CLOSURE_PASS_ZERO_CREDIT_NO_RUNTIME_QUALIFICATION'if report['required_metadata_closure']else'DIAGNOSTIC_UNRESOLVED'
    (x.output/'comparison_PRIVATE.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'paths':len(report['rows']),'current_segment_checks':len(report['current_mapped_segments']),'historical_v2_segments':len(report['historical_v2_segments']),'full_identity_checks_agree':report['all_current_full_identity_checks_agree'],'required_metadata_closure':report['required_metadata_closure']}))
    return 0 if report['required_metadata_closure']else 2

if __name__=='__main__':raise SystemExit(main())
