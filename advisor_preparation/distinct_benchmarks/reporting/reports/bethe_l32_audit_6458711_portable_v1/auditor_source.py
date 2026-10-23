#!/usr/bin/env python3
"""Portable L32 adaptation of the frozen L24 arithmetic auditor; no physics."""
import argparse,hashlib,json,math,os,re,shutil,subprocess,sys,tarfile,tempfile,tomllib
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
checks=0
def require(v,label):
    global checks
    checks+=1
    if not v:raise ValueError(label)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def finite(v):
    if isinstance(v,dict):return all(finite(x) for x in v.values())
    if isinstance(v,list):return all(finite(x) for x in v)
    return not isinstance(v,float) or math.isfinite(v)
def record(p):return {'path':os.path.relpath(p,BASE),'bytes':p.stat().st_size,'sha256':sha(p)}

def main():
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,default=BASE/'results/bethe_l32_pilot_6458711');p.add_argument('--output',type=Path,default=BASE/'reporting/reports/bethe_l32_audit_6458711_portable_v1');a=p.parse_args()
    raw=a.raw.resolve();out=a.output.resolve();source=BASE/'bethe_ansatz_l32';cfgpath=source/'pilot.toml'
    require(not out.exists(),'New immutable output')
    cfg=tomllib.loads(cfgpath.read_text());data=tomllib.loads((raw/'target.toml').read_text());summary=json.loads((raw/'summary.json').read_text());oracle=json.loads((raw/'oracle.json').read_text())
    require(finite(data) and finite(summary) and finite(oracle),'Finite receipts')
    require(cfg['L']==[32] and cfg['phase']=='PILOT_ONLY' and cfg['scientific_completion_allowed'] is False,'L32 pilot only')
    require(data['manifest_sha256']==oracle['manifest_sha256']==summary['manifest_sha256']==sha(cfgpath),'Exact manifest bindings')
    require(data['adapter_sha256']==sha(source/'pilot_target.jl') and oracle['oracle_source_sha256']==sha(source/'oracle_lieb_wu.py'),'Actual adapter/oracle source')
    require(data['status']=='PILOT_ONLY_COMPLETE' and data['scientific_completion_allowed'] is False,'Completed pilot not science')
    for field in ['julia_version'] :require(data[field]==cfg[field],'Native Julia version')
    for field,key in [('ITensors_version','itensors_version'),('ITensorMPS_version','itensormps_version')]:require(data[field]==cfg[key],'Native MPS versions')
    require(oracle['versions']=={k:cfg[k+'_version'] for k in ['numpy','scipy']},'Oracle numerical versions')
    user=BASE.parent/'validated_package/HubbardRenyiDQMC';native=[]
    for cf,tf,rel in [('source_sha256','source_sha256','analysis/compute_s2_itensor_dmrg.jl'),('project_sha256','project_sha256','Project.toml'),('manifest_sha256','environment_manifest_sha256','Manifest.toml')]:
        f=user/rel;require(sha(f)==cfg[cf]==data[tf],'User source/environment '+rel);native.append(f)
    freeze=BASE/'source_releases/bethe-l32-pilot-v3-freeze.json';fr=json.loads(freeze.read_text());archive=freeze.with_name('bethe-l32-pilot-v3-source.tar.gz')
    require(sha(archive)==fr['archive_sha256']=='f89d0b905a57988154d6c53dac6dd5ff307ac4b60848dd89775f53711816c581','Frozen source archive')
    with tarfile.open(archive) as tf:
        members=tf.getmembers();require(len(members)==len(fr['files'])==12,'Full12-file source closure')
        for m in members:
            require(m.isfile() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts,'Safe regular archive member')
            value=tf.extractfile(m).read();require(hashlib.sha256(value).hexdigest()==fr['files'][m.name] and value==(BASE/m.name).read_bytes(),'Archive/local frozen byte identity')
    prior_tool=BASE/'reporting/reports/bethe_l24_audit_6458547_portable_v2/auditor_source.py';prior_report=prior_tool.with_name('audit.json')
    require(sha(prior_tool)==json.loads(prior_report.read_text())['auditor_sha256'],'Preserved portable L24 auditor')
    ep=cfg['engineering_prerequisite'];prior=source/ep['audit_relative'];anchor=source/ep['anchor_relative']
    require(sha(prior)==ep['audit_sha256'] and sha(anchor)==ep['preexecution_engineering_anchor_sha256'],'Preserved prior engineering evidence hashes')
    previous=json.loads(prior.read_text());eligible=previous['next_L32_engineering_eligibility'];anchor_data=eligible['preexecution_anchor']
    require(previous['status']=='PASS' and previous['checks']==168 and previous['scientific_completions']==0,'Prior168-check engineering audit')
    require(anchor_data['commit']==ep['preexecution_engineering_anchor_commit']=='66950c791f7f939ca176e5427be62cd8dfca8348','Exact published engineering anchor commit')
    ab=anchor.read_bytes();require(hashlib.sha1(b'blob '+str(len(ab)).encode()+b'\0'+ab).hexdigest()==anchor_data['git_blob_sha']=='6f09e3c24b16752df0bf6577ae9d9c631dc20af7','Published anchor Git blob')
    require('at most 30 minutes and peak RSS at most 2300 MiB' in ab.decode(),'Prior predeclared admission thresholds')
    require(eligible['job_elapsed_seconds']==ep['L24_job_elapsed_seconds']==991<=ep['maximum_L24_wall_seconds'],'Actual prerequisite elapsed')
    require(eligible['peak_rss_kib']==ep['L24_peak_rss_kib']==1484548<=ep['maximum_L24_peak_rss_MiB']*1024,'Actual prerequisite memory')
    remote_path=raw/'remote_transfer_metadata.json';remote=json.loads(remote_path.read_text());expected={f'results/bethe_l32_pilot_6458711/L32_{s}_{arm}.toml' for s in cfg['sectors'] for arm in cfg['arms']}
    expected|={f'results/bethe_l32_pilot_6458711/{f}' for f in ['target.toml','oracle.json','summary.json','oracle.resources.txt','target.resources.txt']};expected.add('logs/bethe_l32_pilot_6458711.log')
    require(remote['job_id']=='6458711' and set(remote['files'])==expected,'Exactly12 named remote text transfers')
    for rel,v in remote['files'].items():
        f=BASE/rel;require(f.is_file() and not f.is_symlink() and f.stat().st_size==v['bytes'] and sha(f)==v['sha256'],'Remote/local byte and SHA identity '+rel)
    jobs=[s.split('|') for s in remote['sacct'].splitlines() if s.startswith('6458711|')];require(len(jobs)==1,'Unique job scheduler record');j=jobs[0]
    require(j[:6]==['6458711','COMPLETED','0:0','00:21:56','1','3500M'],'Completed scheduler/resources')
    require(j[6]=='00:45:00' and cfg['resources']['wall_cap_minutes']==60,'Preserved administrative wall-limit override60to45')
    require(j[8]=='2026-10-22T16:22:48','Completion before user stop deadline')
    rows={(r['sector'],r['arm']):r for r in data['rows']};require(len(rows)==len(data['rows'])==6 and set(rows)=={(s,a) for s in cfg['sectors'] for a in cfg['arms']},'All six sector/precision solves')
    g=cfg['gates']
    for (sector,arm),r in rows.items():
        require(tomllib.loads((raw/f'L32_{sector}_{arm}.toml').read_text())==r,'Individual solve/aggregate exact identity')
        ac=cfg['arms'][arm];nu,nd=16+int(sector=='plus'),16-int(sector=='minus')
        require((r['L'],r['Nup_expected'],r['Ndown_expected'])==(32,nu,nd),'Expected L32 charge/spin sector')
        require(r['seed']==ac['seed_base']+32000+10*(cfg['sectors'].index(sector)+1),'Unchanged declared seed algorithm')
        require(abs(r['norm']-1)<=g['norm_abs'],'Norm')
        require(abs(r['Nup']-nu)<=g['sector_abs'] and abs(r['Ndown']-nd)<=g['sector_abs'],'Actual particle/spin sector')
        require(abs(r['energy']-r['reported_energy'])<=g['reported_contracted_energy_abs'],'Reported/contracted physical energy')
        require(g['variance_min']<=r['variance'] and r['variance']/32<=ac['variance_per_site_max'],'Unchanged variance gates')
        require(abs(r['variance_per_site']-r['variance']/32)<1e-12,'Variance units')
        require(abs(r['energy_per_site']-r['energy']/32)<1e-12,'Energy units')
        require(len(r['sweep_truncerrors'])==ac['nsweeps'],'All declared sweep errors retained')
        require(r['final_truncerr']==r['sweep_truncerrors'][-1] and 0<=r['final_truncerr']<=ac['final_truncerr_max'],'Unchanged final truncation gates')
        require(0<r['max_bond']<=ac['chi_max'],'Bond cap')
        if sector=='half':require(g['double_occupancy_min']<=r['double_occupancy_per_site']<=g['double_occupancy_max'],'Half-filled double-occupancy range')
        for key,total in [('site_Nup',r['Nup']),('site_Ndown',r['Ndown']),('site_double_occupancy',32*r['double_occupancy_per_site'])]:require(len(r[key])==32 and abs(sum(r[key])-total)<=1e-9,'All32 site values and normalization '+key)
    values={}
    for arm in cfg['arms']:
        h,minus,plus=(rows[s,arm] for s in ['half','minus','plus']);require(abs(plus['energy']-minus['energy']-cfg['U'])<=g['particle_hole_energy_abs'],'Particle-hole charge-sector identity')
        values[arm]={'energy_per_site':h['energy']/32,'charge_gap':plus['energy']+minus['energy']-2*h['energy'],'double_occupancy':sum(h['site_double_occupancy'])/32}
    errors={k:abs(values['coarse'][k]-values['refined'][k]) for k in values['coarse']}
    for key,v in errors.items():require(v<=g['coarse_refined_'+key+'_abs'],'Unchanged refinement stability '+key)
    require(oracle['status']=='COMPLETED' and oracle['oracle_gate']=='PASS','Complete independent thermodynamic reference')
    for name,value in oracle['values'].items():
        require(value['estimated_absolute_error']<=cfg['oracle']['maximum_estimated_absolute_error'],'Oracle error budget '+name)
        require(all(abs(value['adaptive_integral']-v)<=cfg['oracle']['maximum_estimated_absolute_error'] for v in value['fixed_gauss'].values()),'Independent quadrature agreement '+name)
    require(oracle['energy_derivative_error']<=cfg['oracle']['energy_derivative_check_abs'],'Energy derivative/double-occupancy identity')
    with tempfile.TemporaryDirectory(prefix='bethe_l32_postprocess_') as temp:
        tr=Path(temp)
        with tarfile.open(archive) as tf:
            for m in tf:
                dest=tr/m.name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(tf.extractfile(m).read())
        replay=tr/'replay.json';fs=tr/'bethe_ansatz_l32'
        subprocess.run([sys.executable,str(fs/'summarize_pilot.py'),'--manifest',str(fs/'pilot.toml'),'--oracle',str(raw/'oracle.json'),'--target',str(raw/'target.toml'),'--output',str(replay)],check=True,capture_output=True,text=True)
        rb=replay.read_bytes();require(rb==(raw/'summary.json').read_bytes(),'Byte-identical fresh frozen summarizer replay')
    require(summary['pilot_quality']=='PILOT_QUALITY_PASS' and summary['distinct_ideas_passed']==0 and summary['solves']==6,'Pilot quality with zero scientific completions')
    require({r['id'] for r in summary['rows']}==set(cfg['catalogue_ids']) and all(r['status']=='UNRESOLVED' for r in summary['rows']),'Raw summary retains unexecuted thermodynamic scope')
    resources=(raw/'target.resources.txt').read_text();rss=int(re.search(r'Maximum resident set size \(kbytes\): (\d+)',resources).group(1))
    require('Exit status: 0' in resources and rss<=cfg['resources']['memory_mb']*1024,'Successful bounded process memory')
    require(data['elapsed_seconds']<cfg['resources']['target_command_cap_seconds'] and j[3]=='00:21:56','Actual runtime within unchanged command and administrative allocation caps')
    inputs=[BASE/rel for rel in remote['files']]+[remote_path,archive,freeze,prior_tool,prior_report]+[BASE/rel for rel in fr['files']]+native
    records=[record(f) for f in dict.fromkeys(inputs)];out.mkdir();(out/'frozen_summary_replay.json').write_bytes(rb);shutil.copyfile(__file__,out/'auditor_source.py')
    report={'schema':'bethe_l32_portable_pilot_audit_v1','status':'PASS','checks':checks,'scientific_completions':0,'catalogue_status':{i:'NOT_RUN' for i in cfg['catalogue_ids']},'refined':values['refined'],'coarse_refined_errors':errors,'source_archive_sha256':sha(archive),'source_members':12,'remote_text_files':12,'native_source_environment_pins':3,'auditor_sha256':sha(__file__),'adapted_from_portable_L24_auditor_sha256':sha(prior_tool),'frozen_summary_replay_sha256':hashlib.sha256(rb).hexdigest(),'resources':{'job_elapsed_seconds':1316,'target_adapter_seconds':data['elapsed_seconds'],'batch_peak_rss_kib':1735276,'target_process_peak_rss_kib':rss,'frozen_source_wall_minutes':60,'administrative_effective_wall_minutes':45,'administrative_override_changes_physics_or_gates':False},'inputs':records,'input_sha256':{r['path']:r['sha256'] for r in records},'physics_executed':False,'future_growth_planned':False,'scope':'Postprocessing only; single finiteL32 OBC pilot is not a thermodynamic benchmark. A18/A19/A20 remainNOT_RUN; no new scientific jobs or submissions. Prior engineering admission evidence preserved without deriving permission for larger execution.'}
    (out/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    (out/'README.md').write_text('# Portable L32 pilot audit\n\nAll six native DMRG sector/precision solves pass the unchanged pilot-quality gates. Every one of12 remote text/log transfers and12 frozen source members matches its exactSHA256 and byte count. The frozen stdlib summarizer replay is byte-identical. No model, optimizer, eigensolver or Monte Carlo was rerun.\n\nA18/A19/A20 remain NOT_RUN, with0 scientific completions. The32-site finite-size results are not thermodynamic Lieb–Wu equality claims. All previous source, engineering-anchor and L24 audit evidence remains unchanged. No growth planning or new submission is made.\n\nThe frozen source requested60minutes; the effective scheduler allocation was administratively shortened to45minutes under the user stop instruction. Both values remain in the receipt. The job completed at16:22:48EDT in21m56s. Target process peakRSS was'+str(rss)+'KiB; scheduler batch peakRSS was1735276KiB. These are distinct recorded measures.\n')
    files=[f for f in out.iterdir() if f.is_file()];(out/'SHA256SUMS').write_text(''.join(sha(f)+'  '+f.name+'\n' for f in sorted(files)))
    print(json.dumps({'status':'PASS','checks':checks,'scientific_completions':0,'refined':values['refined'],'coarse_refined_errors':errors,'audit_sha256':sha(out/'audit.json')},indent=2))
if __name__=='__main__':main()
