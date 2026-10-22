#!/usr/bin/env python3
"""Independent stdlib audit of immutable A15/A16/A17 receipts; no physics run."""
import argparse, hashlib, itertools, json, math, pathlib, subprocess, sys, tarfile, tempfile, tomllib

BASE = pathlib.Path(__file__).resolve().parents[1]
def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def require(value, message):
    if not value: raise ValueError(message)
def finite(x):
    if isinstance(x, float): require(math.isfinite(x), 'Nonfinite numeric receipt')
    elif isinstance(x, dict):
        for v in x.values(): finite(v)
    elif isinstance(x, list):
        for v in x: finite(v)
def err(a,b):
    if isinstance(a,list) or isinstance(b,list):
        require(isinstance(a,list) and isinstance(b,list) and len(a)==len(b), 'Array shape')
        return max((err(x,y) for x,y in zip(a,b)),default=0.)
    return abs(a-b)
def matrix(flat,n):
    require(len(flat)==n*n,'Full flattened matrix');return [flat[i*n:(i+1)*n] for i in range(n)]
def mul(a,b): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def kron(a,b):return [[a[i//len(b)][j//len(b[0])]*b[i%len(b)][j%len(b[0])] for j in range(len(a[0])*len(b[0]))] for i in range(len(a)*len(b))]
def identity(n):return [[float(i==j) for j in range(n)] for i in range(n)]
def keyed(rows,fields,expected):
    d={tuple(r[k] for k in fields):r for r in rows}
    require(len(d)==len(rows) and set(d)==set(expected),'Complete unique fixture grid');return d
def record(path):return {'path':str(path.relative_to(BASE)) if path.is_relative_to(BASE) else path.name,'bytes':path.stat().st_size,'sha256':sha(path)}

def audit(a):
    folder=a.results.resolve();source=a.source.resolve();freeze=a.freeze.resolve()
    cfg=tomllib.loads((source/'manifest.toml').read_text());o=json.loads((folder/'oracle.json').read_text())
    t=tomllib.loads((folder/'target.toml').read_text());c=json.loads((folder/'comparison.json').read_text())
    for r in (o,t,c):finite(r)
    meta=json.loads((folder/'remote_transfer_metadata.json').read_text())
    require(meta['job_id']=='6458671' and '6458671|COMPLETED|0:0|00:02:05|' in meta['sacct'],'Completed scheduler receipt')
    for rel,want in meta['files'].items():
        path=BASE/rel;require(path.is_file() and not path.is_symlink(),'Named text file')
        require(path.stat().st_size==want['bytes'] and sha(path)==want['sha256'],'Remote/local transfer identity '+rel)
    required={'comparison.json','end_utc.txt','oracle.json','oracle.resources.txt','start_utc.txt','target.resources.txt','target.toml'}
    rows=[x.split(None,1) for x in (folder/'receipts.sha256').read_text().splitlines()]
    require({pathlib.Path(x[1]).name for x in rows}==required and len(rows)==7,'Complete raw checksum receipt')
    for h,name in rows:require(sha(folder/pathlib.Path(name).name)==h,'Job-produced raw checksum')
    fr=json.loads(freeze.read_text());archive=freeze.with_name('two-orbital-analytic-v1-source.tar.gz')
    require(sha(archive)==fr['archive_sha256']=='3ac8111fcbf2b04e8201cd4602da80f55b94b672ff335562942671cb3c06d983','Frozen source archive')
    require(len(fr['files'])==14,'Exact frozen source closure')
    with tarfile.open(archive) as tf:
        members=[m for m in tf.getmembers() if m.isfile()]
        require(len(members)==14 and {m.name for m in members}==set(fr['files']),'Exact14 archive members')
        for member in members:
            require(not member.issym() and not member.islnk(),'No source link')
            raw=tf.extractfile(member).read();local=source.parent/member.name
            require(hashlib.sha256(raw).hexdigest()==fr['files'][member.name] and raw==local.read_bytes(),'Archive/local source bytes')
    provenance=json.loads((source/'source_provenance.json').read_text())
    original=(source/provenance['source_archive']).resolve()
    require(sha(original)==provenance['archive_sha256'],'Original user source archive')
    with tarfile.open(original) as tf:
        for r in provenance['files']:
            raw=tf.extractfile(r['archive_member']).read()
            require(len(raw)==r['bytes'] and hashlib.sha256(raw).hexdigest()==r['sha256']==cfg['source_hashes'][r['path']],'Original user source pin')
            require(raw==(source/cfg['source_root']/r['path']).read_bytes(),'Verbatim user source snapshot')
    require(t['mode']=='full' and t['release']==o['release']==cfg['release'],'Full immutable release')
    require(t['manifest_sha256']==o['manifest_sha256']==c['manifest_sha256']==sha(source/'manifest.toml'),'Manifest identity')
    require(t['driver_sha256']==sha(source/'target.jl') and o['driver_sha256']==sha(source/'oracle.py'),'Actual driver identities')
    require(t['source_hashes']==cfg['source_hashes'],'Actual native user source pins')
    require(t['julia_version']==cfg['julia_version'] and t['blas_threads']==t['julia_threads']==1,'Pinned Julia/thread runtime')
    require(o['numpy_version']==cfg['numpy_version'] and o['scipy_version']==cfg['scipy_version'],'Pinned numerical oracle versions')
    require(tuple(map(int,o['python_version'].split('.')[:2]))>=tuple(map(int,cfg['python_minimum'].split('.'))),'Python runtime')
    require(t['slurm_job_id']=='6458671','Target job attribution')
    require(c['target_sha256']==sha(folder/'target.toml') and c['oracle_sha256']==sha(folder/'oracle.json') and c['comparator_sha256']==sha(source/'compare.py'),'Comparator raw bindings')
    with tempfile.TemporaryDirectory(prefix='two_orbital_postprocess_') as td:
        out=pathlib.Path(td)/'comparison.json'
        subprocess.run([sys.executable,str(source/'compare.py'),'--target',str(folder/'target.toml'),'--oracle',str(folder/'oracle.json'),'--output',str(out)],check=True,capture_output=True,text=True)
        require(out.read_bytes()==(folder/'comparison.json').read_bytes(),'Byte-identical frozen comparator replay')
        require(out.read_bytes()==(folder/'comparison_independent_replay_attempt01.json').read_bytes(),'Preserved replay identity')
    statuses={r['id']:r['status'] for r in c['benchmarks']}
    require(len(c['benchmarks'])==3 and statuses=={'A15':'PASS','A16':'PASS','A17':'PASS'},'Exactly three eligible ideas')
    require(c['unsupported']==[dict(id='A25',status='NOT_RUN',reason=cfg['A25']['reason'])],'A25 remainsNOT_RUN')
    require(all(r['checks'] and all(v['passed'] for v in r['checks']) for r in c['benchmarks']),'All796 frozen scientific and eligibility checks')

    # A15: independently reconstruct the declared Kronecker one-body matrix,
    # then check complete external eigenspace projectors, without diagonalizing.
    p=cfg['A15'];eye=identity(2);x=[[0.,1.],[1.,0.]];j=[[0.,1.],[-1.,0.]]
    ax=kron(eye,x);ay=kron(x,eye);hops=p['hoppings'];t1,t2,t3,t4=hops
    terms=[kron(ax,[[t1,0.],[0.,t2]]),kron(ay,[[t2,0.],[0.,t1]]),kron(mul(ax,ay),eye),kron(kron(j,j),x)]
    H=[[-terms[0][i][k]-terms[1][i][k]-t3*terms[2][i][k]+t4*terms[3][i][k] for k in range(8)] for i in range(8)]
    groups=o['A15']['groups'];require(len(o['A15']['energies'])==8 and sum(g['multiplicity'] for g in groups)==8,'All8 one-electron states')
    total=[[0.]*8 for _ in range(8)];projector_residual=0.
    for g in groups:
        P=matrix(g['projector'],8);energy=sum(g['energies'])/g['multiplicity'];HP=mul(H,P)
        require(len(g['energies'])==g['multiplicity'] and max(g['energies'])-min(g['energies'])<=p['degeneracy_tolerance'],'Whole degeneracy group')
        projector_residual=max(projector_residual,err(HP,[[energy*z for z in row] for row in P]))
        require(err(mul(P,P),P)<=p['projector_tolerance'] and err(P,list(map(list,zip(*P))))<=p['projector_tolerance'],'Hermitian idempotent projector')
        require(abs(sum(P[i][i] for i in range(8))-g['multiplicity'])<=p['projector_tolerance'],'Projector rank')
        require(err([P[i][i] for i in range(8)],g['site_probabilities'])<=p['projector_tolerance'],'Site projector probabilities')
        require(err([sum(P[i][i] for i in range(a,8,2)) for a in range(2)],g['orbital_weights'])<=p['projector_tolerance'],'Orbital projector weights')
        for i in range(8):
            for k in range(8):total[i][k]+=P[i][k]
    require(projector_residual<=p['projector_tolerance'] and err(total,identity(8))<=p['projector_tolerance'],'Complete independent spectral resolution')
    keyed(t['A15'],['U','Nup','Ndown'],[(u,*s) for u in p['interactions'] for s in p['spin_sectors']])

    # A16: recompute areas, peak positions and the observable's exact discrete
    # Gaussian broadening from raw continuum bin masses; no band sampling.
    p=cfg['A16'];require(p['meshes']==[256,512,1024,2048] and p['held_out_mesh']==2048,'Frozen held-out mesh protocol')
    meshes=keyed(t['A16'],['Nk'],[(n,) for n in p['meshes']]);de=(p['energy_max']-p['energy_min'])/p['bins']
    energy=[p['energy_min']+(i+.5)*de for i in range(p['bins'])]
    masses=o['A16']['bin_masses'];require(len(masses)==512 and min(masses)>0 and abs(sum(masses)-1)<=p['normalization_tolerance'],'Positive normalized continuum bin masses')
    oracle_curves=keyed(o['A16']['curves'],['sigma'],[(s,) for s in p['smoothing_widths']]);convolution_error=0.;area_error=0.;peak_error=0.
    for (sigma,),curve in oracle_curves.items():
        radius=min(p['bins']-1,math.ceil(4*sigma/de));kernel=[math.exp(-.5*(k*de/sigma)**2) for k in range(-radius,radius+1)];z=sum(kernel);kernel=[v/z for v in kernel]
        density=[sum(kernel[q+radius]*masses[i-q]/de for q in range(-radius,radius+1) if 0<=i-q<len(masses)) for i in range(len(masses))]
        z=sum(density)*de;density=[v/z for v in density]
        convolution_error=max(convolution_error,err(density,curve['density']));require(convolution_error<1e-11,'Continuum raw-bin Gaussian observable replay')
    for (nk,),row in meshes.items():
        require(row['total_counts']==2*nk*nk and abs(row['bin_width']-de)<1e-14 and err(row['energy'],energy)<1e-14,'Complete native two-band DOS normalization/grid')
        curves=keyed(row['curves'],['sigma'],[(s,) for s in p['smoothing_widths']])
        for curve in curves.values():
            require(len(curve['density'])==512,'All512 DOS bins')
            area=sum(curve['density'])*de;imax=max(range(512),key=curve['density'].__getitem__)
            area_error=max(area_error,abs(area-1),abs(area-curve['normalization']))
            peak_error=max(peak_error,abs(curve['peak']-curve['density'][imax]),abs(curve['peak_energy']-energy[imax]))
    require(area_error<=p['normalization_tolerance'] and peak_error<=1e-12,'Recomputed DOS area and actual peak metadata')
    require(abs(sum(masses[:16])-o['A16']['lower_edge_mass'])<=1e-14 and abs(sum(masses[-16:])-o['A16']['upper_edge_mass'])<=1e-14,'Edge-window bin convention')

    # A17: recalculate all27 product-state atom fixtures directly from the
    # four local Boltzmann weights, not from either producer's observable code.
    p=cfg['A17'];grid=list(itertools.product(p['interactions'],p['centered_chemical_potentials'],p['inverse_temperatures']))
    atoms=keyed(t['A17'],['U','mu','beta'],grid);keyed(o['A17'],['U','mu','beta'],grid);atom_errors={}
    for (u,mu,beta),row in atoms.items():
        energies=[u/4,-u/4-mu,-u/4-mu,u/4-2*mu];low=min(energies);w=[math.exp(-beta*(e-low)) for e in energies];z=sum(w);w=[x/z for x in w]
        mean_n=w[1]+w[2]+2*w[3];moment=(w[1]+w[2])/4
        probabilities=[]
        for state in range(16):
            probabilities.append(w[((state>>0)&1)+2*((state>>2)&1)]*w[((state>>1)&1)+2*((state>>3)&1)])
        analytic={'probabilities':probabilities,'density':2*mean_n,'double_occupancy':2*w[3],
                  'density_connected':0.,'spin_connected':[0.]*9,'local_spin_squared':[moment]*3,
                  'local_density_variance':w[1]+w[2]+4*w[3]-mean_n*mean_n,'S_alm_q0':2*moment,'delta_alm':0.}
        require(row['dimension']==16 and row['sector_count']==9 and len(row['spin_connected'])==9,'Complete16-state/9-sector/all9-spin atom trace')
        for name,want in analytic.items():
            delta=err(row[name],want);atom_errors[name]=max(atom_errors.get(name,0.),delta)
            require(delta<=p['observable_tolerance'],'Independent atomic product formula '+name)
        for name in ['density','double_occupancy','S_alm_q0','delta_alm']:
            delta=err(row['native_'+name],analytic[name]);atom_errors['native_'+name]=max(atom_errors.get('native_'+name,0.),delta)
            require(delta<=p['observable_tolerance'],'Actual native atom helper '+name)
    metrics={}
    for r in c['benchmarks']:
        checks=r['checks'];metrics[r['id']]={'status':'PASS','frozen_check_count':len(checks),'failed_checks':0}
        if r['id']=='A15':
            for suffix in ['physical_energies','projector','site_probabilities','orbital_weights','native_ground_energy']:
                metrics['A15']['max_'+suffix+'_error']=max(x['value'] for x in checks if x['name'].endswith(suffix))
        elif r['id']=='A16':
            metrics['A16']['metric_gates']=[x for x in checks if any(s in x['name'] for s in ['full_curve','last_mesh','edge_mass','peak','oracle_'])]
    metrics['A15'].update(fixture_count=6,independent_H_projector_residual_max=projector_residual)
    metrics['A16'].update(meshes=[256,512,1024,2048],held_out_mesh=2048,curves_per_mesh=2,bins_per_curve=512,independent_raw_mass_convolution_max_error=convolution_error,recomputed_area_max_error=area_error,recomputed_peak_metadata_max_error=peak_error)
    metrics['A17'].update(fixture_count=27,sector_count_per_fixture=9,states_per_fixture=16,spin_connected_components_per_fixture=9,independent_analytic_errors=atom_errors)
    resources={}
    for which in ['target','oracle']:
        text=(folder/(which+'.resources.txt')).read_text();require('Exit status: 0' in text,'Successful timed command')
        resources[which]={'elapsed_wall_text':next(x.split(': ',1)[1] for x in text.splitlines() if 'Elapsed (wall clock)' in x),'maximum_rss_kib':int(next(x.split(': ',1)[1] for x in text.splitlines() if 'Maximum resident set size' in x)),'adapter_compute_elapsed_seconds':t['elapsed_seconds'] if which=='target' else o['elapsed_seconds']}
    resources['scheduler_elapsed_seconds']=125
    input_paths=[folder/n for n in sorted(required|{'receipts.sha256','remote_transfer_metadata.json','comparison_independent_replay_attempt01.json'})]+[BASE/'logs/two_orbital_analytic_full_6458671.log',freeze,archive,original]+[source.parent/rel for rel in fr['files']]
    inputs=[record(p) for p in input_paths]
    return {'status':'PASS','scientific_completions':3,'catalogue_ids':['A15','A16','A17'],'source_archive_sha256':sha(archive),'benchmarks':metrics,'unsupported':[{'id':'A25','status':'NOT_RUN'}],
        'source_closure':{'frozen_members':14,'verbatim_user_source_members':4,'archive_sha256':sha(archive),'original_source_archive_sha256':sha(original)},
        'independent_comparator_replay':'BYTE_IDENTICAL','resources':resources,'inputs':inputs,'input_sha256':{r['path']:r['sha256'] for r in inputs},
        'limits':'Audit replays immutable postprocessing and independently recomputes spectral-projector identities, raw-bin broadening and atomic formulas. It executes no band mesh, eigensolver, Monte Carlo or new model. A16 continuum integral generation remains the pinned independent SciPy oracle, whose mathematical normalization and reduced-BZ identity were reviewed.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--results',type=pathlib.Path,default=BASE/'results/two_orbital_analytic_full_6458671');p.add_argument('--source',type=pathlib.Path,default=BASE/'two_orbital_analytic');p.add_argument('--freeze',type=pathlib.Path,default=BASE/'source_releases/two-orbital-analytic-v1-freeze.json');p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
    require(not a.output.exists(),'Use a new immutable audit output')
    try:result=audit(a)
    except Exception as exc:result={'status':'UNRESOLVED','scientific_completions':0,'reason':f'{type(exc).__name__}: {exc}'}
    result.update(schema='two_orbital_independent_receipt_audit_v1',auditor_sha256=sha(__file__),no_physics_executed=True)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'scientific_completions':result['scientific_completions'],'benchmarks':{k:v['status'] for k,v in result.get('benchmarks',{}).items()},'reason':result.get('reason')},indent=2))
    if result['status']!='PASS':raise SystemExit(1)
if __name__=='__main__':main()
