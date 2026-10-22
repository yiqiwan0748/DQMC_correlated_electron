#!/usr/bin/env python3
"""Audit native ALF fixed-configuration exports. This is not P21 science.

Run inside the pinned QuSpin NumPy/SciPy environment on ISAAC. The oracle in
this adapter proof is explicit dense algebra; the future external reference
for P21 remains native interacting finite-temperature ALF measurements.
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def read_native_bins(path, width):
    rows=[]
    for lineno,line in enumerate(Path(path).read_text().splitlines(),1):
        tokens=re.findall(r'\(([^()]*)\)',line)
        if len(tokens)!=width:raise ValueError(f'{path}:{lineno}: incomplete native vector')
        prefix=line[:line.index('(')].strip()
        if int(prefix)!=width+1:raise ValueError('Native vector width marker differs')
        nums=[complex(*map(float,x.replace('D','E').split(','))) for x in tokens]
        tail=line[line.rfind(')')+1:].strip()
        denominator=float(tail.replace('D','E'))
        rows.append((nums,denominator))
    if not rows:raise ValueError('No complete native bins')
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',type=Path,required=True)
    p.add_argument('--alf-root',type=Path,required=True)
    p.add_argument('--overlay-receipt',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise ValueError('Choose a new proof receipt')
    import numpy as np
    import scipy
    from scipy.linalg import expm
    here=Path(__file__).resolve().parent
    overlay=json.loads(a.overlay_receipt.read_text())
    source={r['path']:r['sha256'] for r in json.loads((here/'upstream_files.json').read_text())['files'] if r['path'].startswith('Prog/')}
    source.update({r['path']:r['observer_sha256'] for r in overlay['files'] if r['origin']=='alf'})
    for path,want in source.items():
        if sha(a.alf_root/path)!=want:raise ValueError('Actual source differs: '+path)
    L=4;M=10;nf=2;dt=.1;U=2.;I=np.eye(L,dtype=complex)
    records={};path=a.run_dir/'P21_native_propagators.tsv'
    for lineno,line in enumerate(path.read_text().splitlines(),1):
        if line.startswith('#') or not line.strip():continue
        kind,t,s,i,j,re,im=line.split()
        key=(kind,int(t),int(s),int(i),int(j))
        if key in records:raise ValueError(f'Duplicate native element at line{lineno}')
        records[key]=complex(float(re.replace('D','E')),float(im.replace('D','E')))
    def matrix(kind,t,s):
        return np.array([[records[kind,t,s,i,j] for j in range(L)] for i in range(L)])
    K=np.zeros((L,L))
    for i in range(L):K[i,(i+1)%L]=-1.;K[(i+1)%L,i]=-1.
    expectedK=expm(-dt*K)
    errs={'HS_operator_identity':0.,'hopping':0.,'coupling':0.,'field_pattern':0.,'VK_slice':0.,'inverse':0.,'G00':0.,'GT0':0.,'G0T':0.,'endpoint':0.}
    alpha=math.acosh(math.exp(dt*U/2.))
    for up in (0,1):
        for down in (0,1):
            q=up-down;centered=math.exp(-dt*U*(up-.5)*(down-.5))
            gaussian=math.exp(dt*U*q*q/2.-dt*U/4.)
            hirsch=math.exp(-dt*U/4.)*math.cosh(alpha*q)
            errs['HS_operator_identity']=max(errs['HS_operator_identity'],abs(gaussian-centered),abs(hirsch-centered))
    wrong_order_discrepancy=0.;Bs={}
    for s in range(nf):
        Bs[s]=[]
        for t in range(1,M+1):
            fields=np.array([records['field',t,s,i,0] for i in range(L)])
            g=np.array([records['coupling',t,s,i,0] for i in range(L)])
            expected_fields=np.array([((3*(i+1)+2*t)%7-3)/5. for i in range(L)])
            errs['field_pattern']=max(errs['field_pattern'],float(np.max(np.abs(fields-expected_fields))))
            errs['coupling']=max(errs['coupling'],float(np.max(np.abs(g-((-1.)**s)*np.sqrt(dt*U)))))
            actualK=matrix('Kexp',t,s);actualB=matrix('B',t,s);actualInv=matrix('Binv',t,s)
            D=np.diag(np.exp(g*fields));expectedB=D@expectedK
            errs['hopping']=max(errs['hopping'],float(np.max(np.abs(actualK-expectedK))))
            errs['VK_slice']=max(errs['VK_slice'],float(np.max(np.abs(actualB-expectedB))))
            errs['inverse']=max(errs['inverse'],float(np.max(np.abs(actualB@actualInv-I))))
            wrong_order_discrepancy=max(wrong_order_discrepancy,float(np.max(np.abs(actualB-expectedK@D))))
            Bs[s].append(actualB)
    n=(M+1)*nf*L*L
    gtrows=read_native_bins(a.run_dir/'P21_GT0_scal',n)
    gzrows=read_native_bins(a.run_dir/'P21_G0T_scal',n)
    counts=read_native_bins(a.run_dir/'P21_Counts_scal',M+1)
    if len(gtrows)!=1 or len(gzrows)!=1 or len(counts)!=1:raise ValueError('Proof requires exactly one native grid')
    if counts[0][1]!=1. or any(x!=1+0j for x in counts[0][0]):raise ValueError('Incomplete/repeated native tau grid')
    if abs(gtrows[0][1])<.999999999999 or abs(gtrows[0][1])>1.000000000001:raise ValueError('Unexpected real proof sign')
    if gtrows[0][1]!=gzrows[0][1]:raise ValueError('Unpaired native denominators')
    gt=np.array(gtrows[0][0]).reshape(M+1,nf,L,L)/gtrows[0][1]
    gz=np.array(gzrows[0][0]).reshape(M+1,nf,L,L)/gzrows[0][1]
    for s in range(nf):
        total=I.copy()
        for b in Bs[s]:total=b@total
        g00=np.linalg.inv(I+total)
        errs['G00']=max(errs['G00'],float(np.max(np.abs(matrix('G00',0,s)-g00))))
        prefix=I.copy()
        for t in range(M+1):
            if t:prefix=Bs[s][t-1]@prefix
            expectedGT=prefix@g00;expectedGZ=(g00-I)@np.linalg.inv(prefix)
            errs['GT0']=max(errs['GT0'],float(np.max(np.abs(gt[t,s]-expectedGT))))
            errs['G0T']=max(errs['G0T'],float(np.max(np.abs(gz[t,s]-expectedGZ))))
        errs['endpoint']=max(errs['endpoint'],float(np.max(np.abs(gt[0,s]-g00))),float(np.max(np.abs(gz[0,s]-(g00-I)))),float(np.max(np.abs(gt[M,s]-(I-g00)))),float(np.max(np.abs(gz[M,s]+g00))))
    expected_count=M*nf*(3*L*L+2*L)+nf*L*L
    if len(records)!=expected_count:raise ValueError('Unexpected/missing native records')
    tolerances={k:(5e-10 if k in ('G00','GT0','G0T','endpoint') else 5e-12) for k in errs}
    passed=all(np.isfinite(v) and v<=tolerances[k] for k,v in errs.items()) and wrong_order_discrepancy>1e-5
    receipt={'schema':'p21_alf_native_fixed_configuration_proof_v1','catalogue_id':'P21',
        'status':'ADAPTER_PROOF_PASS_NOT_EXTERNAL_BENCHMARK' if passed else 'ADAPTER_PROOF_FAIL',
        'alf_commit':overlay['alf_commit'],'source_sha256':source,'binary_sha256':sha(a.binary),
        'observer_receipt_sha256':sha(a.overlay_receipt),'driver_sha256':sha(__file__),
        'numpy_version':np.__version__,'scipy_version':scipy.__version__,
        'max_errors':errs,'tolerances':tolerances,'KV_wrong_order_discrepancy':wrong_order_discrepancy,
        'time_indices':list(range(M+1)),'fixed_configuration':True,'monte_carlo_sweeps':0,
        'inputs_sha256':{p.name:sha(p) for p in [path,a.run_dir/'P21_GT0_scal',a.run_dir/'P21_G0T_scal',a.run_dir/'P21_Counts_scal']}}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
    if not passed:raise SystemExit(1)
if __name__=='__main__':main()
