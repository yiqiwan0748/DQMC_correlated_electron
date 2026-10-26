#!/usr/bin/env python3
"""Local saved-data replay, not a pinned ISAAC runtime replay or new benchmark."""
import ast,datetime,hashlib,importlib.util,json,math,pathlib,platform,time,tomllib
from fractions import Fraction
ROOT=pathlib.Path(__file__).resolve().parents[2]
SNAP=ROOT/'reconciliation_20261025/weekend_state/snapshot'
SOURCE=SNAP/'source/x14_full_independent_saved_audit_v1_pc1_20261024'
ORIGINAL=SNAP/'source/x14_saved_compare_recovery_v1_pc1_20261024/original_source'
RAW=ROOT/'results/x14_full_v3_cancelled_6467601/results'
OUT=pathlib.Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def defs(p):return {n.name:n for n in ast.parse(p.read_text()).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
def eq(a,b):return ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False)
start=time.monotonic();checks=[]
def need(name,cond):
 checks.append(dict(name=name,passed=bool(cond)))
 if not cond:raise ValueError(name)
newdefs=defs(SOURCE/'audit_full.py');olddefs=defs(SOURCE/'provenance/original_pilot_auditor.py')
for n in ('sha','read','canonical','Audit','ratio_term','exact_sum','bounds','value','physical'):need('unchanged independent pilot math helper '+n,eq(newdefs[n],olddefs[n]))
for n in ('__init__','take','sign','root','census'):
 a={n.name:n for n in newdefs['Stream'].body if isinstance(n,ast.FunctionDef)};b={n.name:n for n in olddefs['Stream'].body if isinstance(n,ast.FunctionDef)};need('unchanged independent pilot census helper '+n,eq(a[n],b[n]))
need('physical manifest byte exact first version',sha(ORIGINAL/'manifest.toml')==sha(ORIGINAL/'original_v1/manifest.toml')=='953b483b8550012fc81d1db8b843ec8d867c0a615d0ffa8503add927d93ae1fc')
spec=importlib.util.spec_from_file_location('preserved_independent_saved_math',SOURCE/'audit_full.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
c=tomllib.loads((ORIGINAL/'manifest.toml').read_text());n=load(RAW/'native.json');r=load(RAW/'reference.json');receipt=load(SNAP/'results/x14_full_independent_saved_audit_v1_pc1_20261024/audit.json');comp=load(SNAP/'results/x14_saved_compare_recovery_v1_pc1_20261024/comparison.json')
need('full grid140 rows both inputs',len(n['rows'])==len(r['rows'])==140)
for backend,doc in [('native',n),('quspin',r)]:
 spectra=m.expected_spectra(ORIGINAL,backend,c);wanted=[(s,h) for s in spectra for h in c['etas']];need('complete140 distinct prescribed rows '+backend,len(wanted)==140 and len({(v['id'],v['eta'])for v in doc['rows']})==140)
 for (s,h),row in zip(wanted,doc['rows']):need('full source identity '+backend+'/'+row['id']+'/'+str(h),(row['id'],row['role'],row['branch'],row['groups'],row['eta'])==(s['id'],s['role'],s['branch'],s['groups'],h))
a=m.Audit();proof=m.Stream(RAW/'reference_sign_proof.jsonl',a,c,m.OLD_FREEZE)
for row in r['rows']:proof.census(row)
stats=proof.finish(r['complete_sign_proof']);need('replayed exact point/census statistics',stats==receipt['proof_statistics'])
# Additional independent Fraction construction for EVERY distinct point, separate
# from producer's event rational and the independent LCM verifier construction.
rows_by_fp={row['proof_curve_fingerprint']:row for row in r['rows']};fraction_count=0;maxnumerator=0
for (fp,xhex),(sign,num,den,bits) in proof.cache.items():
 row=rows_by_fp[fp];x=Fraction(float.fromhex(xhex));h=Fraction(row['eta']);value=sum((Fraction(w)*(x-Fraction(e))/((x-Fraction(e))**2+h*h)**2 for e,w in zip(row['poles'],row['weights'])),Fraction())
 expected=-1 if value>0 else 1 if value<0 else 0
 need('third exact rational construction '+str(fraction_count),expected==sign and value==Fraction(num,den));fraction_count+=1;maxnumerator=max(maxnumerator,abs(value.numerator).bit_length(),value.denominator.bit_length())
rows,largest=m.physical(n['rows'],r['rows'],c,a);witnesses=m.full_additional_gates(n['rows'],r['rows'],c,a)
need('full140 metric summaries same saved independent receipt',rows==receipt['row_metrics']);need('entire fullcurve result exact same receipt',largest==receipt['max_independent_full_curve_abs']);need('all8 mandatory witnesses unchanged',witnesses==receipt['mandatory_full_witnesses']==[{k:v for k,v in w.items() if k!='visible_pole_threshold'}for w in comp['physical_shift_witnesses']]);need('all scientific checks pass',all(v['passed']for v in a.checks));need('saved independent proof report full PASS',receipt['status']=='PASS'and receipt['scientific_status']=='PASS' and receipt['independent_mathematics']=='PASS' and receipt['check_count']==len(receipt['checks'])==260842 and all(v['passed'] for v in receipt['checks']))
report=dict(schema='WEEKEND_X14_LOCAL_INDEPENDENT_SAVED_MATH_REPLAY_V1',status='PASS',scope='Local postprocessing of unchanged saved140rows only. Uses preserved independently reviewed pilot census/math with documented full-grid extension; adds third Fraction construction for all18267distinct points. No numerical fixture, producer, ISAAC job, pinned-runtime claim or new benchmark count.',python=platform.python_version(),stdlib_only=True,original_pinned_runtime_replay=False,source_sha256=sha(SOURCE/'audit_full.py'),runner_sha256=sha(pathlib.Path(__file__)),input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [RAW/'native.json',RAW/'reference.json',RAW/'reference_sign_proof.jsonl',ORIGINAL/'manifest.toml',SNAP/'results/x14_full_independent_saved_audit_v1_pc1_20261024/audit.json',SNAP/'results/x14_saved_compare_recovery_v1_pc1_20261024/comparison.json']},checks=checks,check_count=len(checks),preserved_auditor_math_check_count=len(a.checks),proof_statistics=stats,additional_Fraction_points=fraction_count,additional_Fraction_maxbits=maxnumerator,full_curve_abs=largest,witnesses=witnesses,elapsed_seconds=time.monotonic()-start,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(OUT/'x14_local_saved_math_replay.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in report.items()if k in ['status','check_count','preserved_auditor_math_check_count','additional_Fraction_points','elapsed_seconds']},indent=2))
