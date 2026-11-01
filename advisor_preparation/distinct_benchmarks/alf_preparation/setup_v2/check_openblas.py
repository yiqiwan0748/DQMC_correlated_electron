#!/usr/bin/env python3
"""Read native OpenBLAS build metadata; no third-party Python imports or physics."""
import ctypes
import hashlib
import json
import os
from pathlib import Path
import platform
import sys

path=Path(sys.argv[1]).resolve(strict=True)
output=Path(sys.argv[2]);assert not output.exists()
lib=ctypes.CDLL(str(path))
config=lib.openblas_get_config;config.restype=ctypes.c_char_p;config.argtypes=[]
text=config().decode('ascii')
assert text.startswith('OpenBLAS 0.3.12'),text
assert 'USE64BITINT' not in text,'OpenBLAS ILP64 is incompatible with default ALF integers'
assert 'NO_LAPACK' not in text,'This OpenBLAS build excludes LAPACK'
required=['ddot_','dgemm_','zgemm_','zheev_']
for name in required:getattr(lib,name)
threads=lib.openblas_get_num_threads;threads.restype=ctypes.c_int;threads.argtypes=[]
count=threads();assert count==1,('Expected one BLAS thread',count)
result={'status':'OPENBLAS_LP64_METADATA_PASS','provider':str(path),'provider_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'config':text,'threads':count,'required_symbols':required,'python':platform.python_version(),'thread_environment':{k:os.getenv(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']},'scientific_samples':0}
with output.open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))
