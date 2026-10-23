#!/usr/bin/env python3
"""Resource-only pilot assessment. It never promotes a catalogue status."""
import argparse,hashlib,json,pathlib,tomllib
HERE=pathlib.Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--target',required=True);p.add_argument('--max-rss-kib',type=int,required=True);p.add_argument('--output',required=True);a=p.parse_args()
    output=pathlib.Path(a.output)
    if output.exists():raise FileExistsError(output)
    m=tomllib.loads((HERE/'manifest.toml').read_text());t=tomllib.loads(pathlib.Path(a.target).read_text());r=m['resources']
    if t.get('mode')!='pilot' or t.get('scientific_completions')!=0:raise ValueError('Not the zero-completion pilot')
    if t['manifest_sha256']!=sha(HERE/'manifest.toml') or t['driver_sha256']!=sha(HERE/'target.jl'):raise ValueError('Unpinned pilot')
    if t['source_hashes']!=m['source_hashes'] or t['julia_version']!=m['julia_version'] or t['blas_threads']!=1 or t['julia_threads']!=1:raise ValueError('Pilot runtime/source mismatch')
    rows=t['A16']
    if [x['Nk'] for x in rows]!=r['pilot_meshes'] or len(t['A15'])!=6 or len(t['A17'])!=27:raise ValueError('Incomplete pilot grid')
    warm=rows[-1]['elapsed_seconds']
    if not (0<warm<float('inf')) or a.max_rss_kib<=0:raise ValueError('Missing or invalid timing/RSS')
    # Empirical resource projection, not a convergence or rigorous runtime guarantee.
    projection=r['projection_safety_factor']*warm*sum(n*n for n in m['A16']['meshes'])/rows[-1]['Nk']**2+r['projection_startup_reserve_seconds']
    eligible=projection<=r['projected_target_limit_seconds'] and a.max_rss_kib<=r['pilot_rss_limit_kib']
    result=dict(release=m['release'],assessment='RESOURCE_ELIGIBLE_FOR_SEPARATE_REVIEW' if eligible else 'RESOURCE_UNRESOLVED',
                scientific_completions=0,catalogue_statuses={i:'NOT_RUN' for i in ['A15','A16','A17','A25']},
                target_sha256=sha(a.target),manifest_sha256=sha(HERE/'manifest.toml'),
                pilot_512_mesh_seconds=warm,measured_pilot_max_rss_kib=a.max_rss_kib,
                projected_full_target_seconds=projection,projection_is_empirical=True,
                limits=dict(seconds=r['projected_target_limit_seconds'],pilot_rss_kib=r['pilot_rss_limit_kib']),
                limitation='Separate root review/submission required. No full2048 mesh, independent oracle, scientific PASS or guaranteed full memory/runtime follows from this pilot.')
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
