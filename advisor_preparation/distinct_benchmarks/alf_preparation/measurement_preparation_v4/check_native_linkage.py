#!/usr/bin/env python3
"""Require the ALF executable's BLAS/LAPACK SONAMEs to resolve to one pinned provider."""
import hashlib,json,subprocess,sys
from pathlib import Path
binary,provider,output=map(Path,sys.argv[1:])
assert not output.exists();provider=provider.resolve(strict=True)
run=subprocess.run(['ldd',str(binary)],capture_output=True,text=True,timeout=20,check=True)
assert 'not found' not in run.stdout,run.stdout
libraries=[]
for line in run.stdout.splitlines():
    if any(name in line for name in ('libopenblas','libblas','liblapack')):
        assert '=>' in line,line
        path=Path(line.split('=>',1)[1].strip().split()[0]).resolve(strict=True)
        libraries.append({'ldd_line':line.strip(),'realpath':str(path)})
assert libraries and all(Path(row['realpath'])==provider for row in libraries),(libraries,str(provider))
result={'schema':'p21_alf_linkage_audit_v1','status':'PINNED_SINGLE_PROVIDER_PASS',
        'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
        'provider':str(provider),'provider_sha256':hashlib.sha256(provider.read_bytes()).hexdigest(),
        'resolved_libraries':libraries,'ldd':run.stdout,'scientific_samples':0}
output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
