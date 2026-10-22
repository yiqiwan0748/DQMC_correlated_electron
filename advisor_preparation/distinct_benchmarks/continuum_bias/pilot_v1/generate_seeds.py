#!/usr/bin/env python3
"""Static integer-only seed generation; never initializes a scientific backend."""
import hashlib,json,tomllib
from pathlib import Path

MASK=(1<<64)-1
def derived_seed(root,k,rep):
    x=root ^ (((k+1)*0x9e3779b97f4a7c15)&MASK) ^ ((rep*0xbf58476d1ce4e5b9)&MASK)
    z=(x+0x9e3779b97f4a7c15)&MASK
    z=((z^(z>>30))*0xbf58476d1ce4e5b9)&MASK
    z=((z^(z>>27))*0x94d049bb133111eb)&MASK
    return (z^(z>>31))%((1<<63)-1)

def contents(cfg):
    namespace=cfg['seed_namespace'];lines=['namespace = '+json.dumps(namespace),'algorithm = "SHA256 first8bytes big-endian masked to63bits, then unchanged native UInt64 derived_seed for replicas"','']
    actual=[];labels=[]
    for index,dt in enumerate(cfg['pilot_dtau'],1):
        for bundle in range(1,cfg['bundles_per_dtau']+1):
            def seed(stream):
                label=f'{namespace}|dtau_index={index}|bundle={bundle}|stream={stream}'
                value=int.from_bytes(hashlib.sha256(label.encode()).digest()[:8],'big')&((1<<63)-1)
                assert value>0;labels.append(value);return value
            energy,root=seed('energy'),seed('entropy')
            replicas=[[derived_seed(root,k,r) for r in (1,2)] for k in cfg['lambda_indices']]
            actual.extend([energy]+[x for pair in replicas for x in pair])
            lines.extend(['[[bundles]]',f'dtau_index = {index}',f'dtau = {dt!r}',f'bundle = {bundle}',f'energy_seed = {energy}',f'entropy_root_seed = {root}','node_replica_seeds = '+json.dumps(replicas),''])
    assert len(actual)==len(set(actual))==36 and min(actual)>0 and len(labels)==len(set(labels))==8
    return '\n'.join(lines)

if __name__=='__main__':
    here=Path(__file__).resolve().parent;cfg=tomllib.loads((here/'manifest.toml').read_text())
    output=here/cfg['seed_file'];payload=contents(cfg)
    if output.exists():assert output.read_text()==payload,'Frozen seed content differs'
    else:output.write_text(payload)
    print('Static seed closure verified:36 disjoint actual initialization seeds;8 disjoint stream labels.')
