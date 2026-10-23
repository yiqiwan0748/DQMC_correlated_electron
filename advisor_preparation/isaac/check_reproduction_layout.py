#!/usr/bin/env python3
"""Check published reproduction inputs, archive closure and links; starts no jobs."""
from __future__ import annotations
import hashlib,json,re,subprocess,tarfile
from pathlib import Path
from urllib.parse import unquote
ROOT=Path(__file__).resolve().parent.parent
ARCHIVES={
 'baseline-source.tar.gz':'11c3c33adaa70a78b8ab10437fba04651e3951c040bb289c4a3ce7ec293420aa',
 'baseline-consumer-20261022.tar.gz':'d25375c07acd2aea687d38da84dd506490135e2bfc4d63388a13f7fd1928888f',
 'corrected-source.tar.gz':'ae174444854849faf5ddb604a2b4a3e8dc88ccdb70f4fb0d9b9ff065dc9636e6',
}
def digest(b):return hashlib.sha256(b).hexdigest()
def members(file):
 with tarfile.open(file) as t:
  out={}
  for e in t.getmembers():
   p=Path(e.name)
   assert not p.is_absolute() and '..' not in p.parts,e.name
   assert e.isfile() or e.isdir(),e.name
   if e.isfile():out[e.name]=t.extractfile(e).read()
  return out

def main():
 archives=[];blobs={}
 for name,expected in ARCHIVES.items():
  path=ROOT/'validation_source'/name;b=path.read_bytes();assert digest(b)==expected,name
  m=members(path);blobs[name]=m
  archives.append({'path':str(path.relative_to(ROOT)),'sha256':digest(b),'bytes':len(b),'regular_files':len(m),'top_level':sorted({p.split('/')[0] for p in m})})
 base=blobs['baseline-source.tar.gz'];checks=base['SHA256SUMS'].decode().splitlines()
 for line in checks:
  sha,path=line.split('  ',1);assert digest(base[path])==sha,path
 original=blobs['baseline-consumer-20261022.tar.gz'];orig=original['validation_design/run_external_validation.jl']
 corrected=blobs['corrected-source.tar.gz'];corr=corrected['validation_design_corrected/run_external_validation.jl']
 assert 't7critical(.01/100)' in orig.decode() and 'split_rhat' not in orig.decode()
 assert 't7critical(.01/400)' in corr.decode() and 'split_rhat' in corr.decode()
 mismatches=[];ncompared=0
 for name,data in corrected.items():
  match=next(((a,b) for a,b in [('renyi_corrected/','validated_package/HubbardRenyiDQMC/'),('validation_design_corrected/','validation_design/')] if name.startswith(a)),None)
  if match:
   a,b=match;local=ROOT/(b+name[len(a):]);ncompared+=1
   if not local.exists() or local.read_bytes()!=data:mismatches.append({'archive_member':name,'current_local_path':str(local.relative_to(ROOT)),'archive_sha256':digest(data),'current_sha256':digest(local.read_bytes()) if local.exists() else None})
 allowed={'renyi_corrected/analysis/atomic_odd_sweep_proof.jl'}
 assert all(e['archive_member'] in allowed for e in mismatches),mismatches
 overlays=[]
 for rel in ['external_oracles/quspin_export.py','external_oracles/cases.json','isaac/diagnostics/diagnose_atomic_odd_sweep_v2.jl','isaac/diagnostics/benchmark_matched_lazy_grover.jl','isaac/atomic_corrected_resume.sbatch','isaac/timing_matched.sbatch','isaac/logs/quspin_requirements.txt']:
  p=ROOT/rel;b=p.read_bytes();overlays.append({'path':rel,'sha256':digest(b),'bytes':len(b)})
 assert next(e['sha256'] for e in overlays if e['path'].endswith('diagnose_atomic_odd_sweep_v2.jl'))=='402f3614afd71ddb87f0977fcc778d9a314b40ff4bd7828443f28597d7f04dad'
 assert next(e['sha256'] for e in overlays if e['path'].endswith('benchmark_matched_lazy_grover.jl'))=='56ce5e3c395df4aeb1fe0b7654e315a81e280c74b2f65bb0b0ab65a010d06c49'
 shellfiles=sorted((ROOT/'isaac').glob('*.sbatch'))+sorted((ROOT/'isaac').glob('*.sh'))
 for p in shellfiles:subprocess.run(['bash','-n',str(p)],check=True,capture_output=True)
 docs=[ROOT/'isaac/ISAAC_REPRODUCE.md',ROOT/'external_oracles/ATOMIC_SAMPLER_DIAGNOSIS.md',ROOT/'external_oracles/README.md',ROOT/'research_record/README.md',ROOT/'research_record/EVIDENCE_LEDGER.md']
 links=[];missing=[]
 for file in docs:
  for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',file.read_text()):
   target=target.strip('<>')
   if re.match(r'^[A-Za-z][A-Za-z0-9+.-]*:',target) or target.startswith('#'):continue
   path=unquote(target.split('#',1)[0]);resolved=file.parent/path
   item={'document':str(file.relative_to(ROOT)),'target':target,'exists':resolved.exists()};links.append(item)
   if not resolved.exists():missing.append(item)
 assert not missing,missing
 result={'status':'PASS_PUBLISHED_INPUTS_AND_REPRODUCTION_LAYOUT','simulations_or_jobs_started':0,'archives':archives,'baseline_inner_files_verified':len(checks),'original_consumer_sha256':digest(orig),'corrected_consumer_sha256':digest(corr),'consumer_protocols_distinguished':True,'corrected_archive_members_compared':ncompared,'later_local_helper_differences':mismatches,'difference_scope':'Only a later module/world-age loading change in an analysis helper. The frozen source archive is preserved; the guide explicitly stages the versioned v2 ISAAC diagnostic. Physical kernels match the frozen corrected archive.','later_support_overlays':overlays,'bash_syntax_files_checked':len(shellfiles),'relative_links_checked':len(links),'missing_relative_links':missing,'links':links,'legacy_alias_note':'Existing validation_results links are compatibility aliases managed by publication. New documents use canonical isaac/results paths.'}
 result['analysis_script_sha256']=digest(Path(__file__).read_bytes())
 (ROOT/'isaac/reproduction_layout_audit.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['links','later_support_overlays','archives']},indent=2))
if __name__=='__main__':main()
