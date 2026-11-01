#!/usr/bin/env python3
"""Independent X06 receipt audit, using Python's standard library only.

Does not import the production comparator, NumPy or SciPy. Jackknife estimates
are means/standard errors of pseudovalues. Binomial tests use direct log-PMF
sums; exact coverage limits invert binomial tails. Student-t quantiles invert
the incomplete-beta CDF using a continued fraction. No simulation is performed.
"""
import argparse,hashlib,json,math,statistics,tomllib
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def avg(x):return math.fsum(x)/len(x)
def variance(x):
    m=avg(x);return math.fsum((v-m)**2 for v in x)/(len(x)-1)
def jackknife(numerators,denominators):
    B=len(numerators);N=math.fsum(numerators);D=math.fsum(denominators)
    theta=N/D
    pseudo=[B*theta-(B-1)*(N-n)/(D-d) for n,d in zip(numerators,denominators)]
    return avg(pseudo),math.sqrt(variance(pseudo)/B),pseudo
def neff(nums,dens,center):
    x=[a-center*b for a,b in zip(nums,dens)];n=len(x);m=avg(x);y=[v-m for v in x]
    g0=math.fsum(v*v for v in y)/n
    if g0<=math.ulp(1.)*max(1.,max(v*v for v in x)):return float(n)
    correlations=[]
    for lag in range(1,min(n//4,16)+1):
        corr=math.fsum(y[i]*y[i+lag] for i in range(n-lag))/((n-lag)*g0)
        if corr<=0.:break
        correlations.append(corr)
    return n/(1.+2.*math.fsum(correlations))
def split_rhat(chains):
    h=len(chains[0])//2
    xs=[c[:h] for c in chains]+[c[h:] for c in chains]
    W=avg([variance(c) for c in xs]);B=h*variance([avg(c) for c in xs])
    tiny=1e-24*max(1.,max(v*v for c in xs for v in c))
    if W<=tiny:return 1. if B<=tiny else math.inf
    return math.sqrt(max(1.,((h-1)*W+B)/(h*W)))

def beta_cf(a,b,x):
    # Lentz continued fraction for the regularized incomplete beta.
    tiny=1e-300;qab=a+b;qap=a+1.;qam=a-1.
    c=1.;d=1.-qab*x/qap
    if abs(d)<tiny:d=tiny
    d=1./d;h=d
    for m in range(1,10001):
        twice=2*m;aa=m*(b-m)*x/((qam+twice)*(a+twice))
        d=1.+aa*d;c=1.+aa/c
        if abs(d)<tiny:d=tiny
        if abs(c)<tiny:c=tiny
        d=1./d;h*=d*c
        aa=-(a+m)*(qab+m)*x/((a+twice)*(qap+twice))
        d=1.+aa*d;c=1.+aa/c
        if abs(d)<tiny:d=tiny
        if abs(c)<tiny:c=tiny
        d=1./d;delta=d*c;h*=delta
        if abs(delta-1.)<2e-14:return h
    raise ArithmeticError('Incomplete beta failed to converge')
def ibeta(a,b,x):
    if x<=0.:return 0.
    if x>=1.:return 1.
    bt=math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log1p(-x))
    if x<(a+1.)/(a+b+2.):return bt*beta_cf(a,b,x)/a
    return 1.-bt*beta_cf(b,a,1.-x)/b
def student_cdf(x,df):
    tail=.5*ibeta(df/2.,.5,df/(df+x*x))
    return 1.-tail if x>=0 else tail
def student_quantile(q,df):
    lo,hi=0.,16.
    for _ in range(100):
        mid=(lo+hi)/2.
        if student_cdf(mid,df)<q:lo=mid
        else:hi=mid
    return (lo+hi)/2.
def logpmf(k,n,p):
    if p==0.:return 0. if k==0 else -math.inf
    if p==1.:return 0. if k==n else -math.inf
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)+k*math.log(p)+(n-k)*math.log1p(-p)
def binomial_tail(first,last,n,p):
    return math.fsum(math.exp(logpmf(k,n,p)) for k in range(first,last+1))
def exact_pvalue(k,n,p):
    cutoff=logpmf(k,n,p)+math.log1p(1e-7)
    return min(1.,math.fsum(math.exp(logpmf(i,n,p)) for i in range(n+1) if logpmf(i,n,p)<=cutoff))
def clopper_pearson(k,n,confidence):
    tail=(1.-confidence)/2.;lo,hi=0.,1.
    for _ in range(100):
        mid=(lo+hi)/2.
        if binomial_tail(k,n,n,mid)<tail:lo=mid
        else:hi=mid
    lower=0. if k==0 else (lo+hi)/2.
    lo,hi=0.,1.
    for _ in range(100):
        mid=(lo+hi)/2.
        if binomial_tail(0,k,n,mid)>tail:lo=mid
        else:hi=mid
    upper=1. if k==n else (lo+hi)/2.
    return [lower,upper]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',type=Path,required=True)
    p.add_argument('--recorded-summary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();own=Path(__file__).resolve().parents[1]
    if a.output.exists():raise ValueError('Choose a new independent audit path')
    cfg=tomllib.loads((own/'manifest.toml').read_text());data=tomllib.loads(a.raw.read_text())
    recorded=json.loads(a.recorded_summary.read_text());seeds=tomllib.loads((own/'seeds.toml').read_text())
    frozen=json.loads((own/'frozen_inputs.json').read_text())
    for f,digest in frozen['files'].items():assert sha(own/f)==digest,(f,'freeze mismatch')
    meta=data['metadata'];ref=cfg['oracle_reference'];N=cfg['final_replicates'];C=cfg['chains_per_replicate'];B=cfg['bins']
    assert meta['catalogue_id']=='X06' and meta['mode']=='final'
    assert meta['manifest_sha256']==sha(own/'manifest.toml') and meta['seeds_sha256']==sha(own/'seeds.toml')
    assert meta['driver_sha256']==sha(own/'run_coverage.jl') and meta['source_hashes']==tomllib.loads((own/'source_hashes.toml').read_text())['files']
    assert meta['reference']==ref and meta['oracle_sha256']==cfg['oracle_sha256']
    assert meta['julia_version']=='1.12.1' and meta['blas_threads']==1
    assert data['completion']['status']=='SHARD_COMPLETE' and data['completion']['replicate_count']==N
    assert len(data['replicates'])==N and [r['id'] for r in data['replicates']]==list(range(1,N+1))
    critical=student_quantile(.975,C-1);assert abs(critical-cfg['student_t7_critical'])<1e-12
    max_replay_error={};allseeds=set();rows=[];all_positive_sign=True
    def compare(key,x,y,tolerance=2e-10):
        error=abs(x-y);max_replay_error[key]=max(max_replay_error.get(key,0.),error)
        assert error<=tolerance*max(1.,abs(x),abs(y)),(key,x,y,error)
    for r in data['replicates']:
        assert r['mode']=='final' and len(r['raw_bins'])==C
        means=[];blockses=[];double=[];eff=[];ps=[]
        for j,raw in enumerate(r['raw_bins'],1):
            seed=int.from_bytes(hashlib.sha256(f"{cfg['seed_namespace']}|mode=final|replicate={r['id']}|chain={j}".encode()).digest()[:8],'big')&((1<<63)-1)
            assert raw['chain']==j and raw['seed']==seed==seeds['final'][str(r['id'])][j-1]
            assert seed not in allseeds;allseeds.add(seed)
            ns,ds=raw['numerator_bins'],raw['denominator_bins']
            assert len(ns)==len(ds)==B and all(math.isfinite(x) for x in ns+ds)
            all_positive_sign&=all(d==cfg['measured_sweeps']/B for d in ds)
            v,se,pseudo=jackknife(ns,ds)
            ns2=[math.fsum(ns[k:k+2]) for k in range(0,B,2)];ds2=[math.fsum(ds[k:k+2]) for k in range(0,B,2)]
            se2=jackknife(ns2,ds2)[1];effective=neff(ns,ds,v)
            means.append(v);blockses.append(se);double.append(se2);eff.append(effective);ps.append(pseudo)
            for key,x,y in [('chain_mean',v,r['chain_means'][j-1]),('chain_block_se',se,r['chain_block_se'][j-1]),('chain_double_se',se2,r['chain_double_block_se'][j-1]),('chain_neff',effective,r['chain_effective_bins'][j-1])]:compare(key,x,y)
        value=avg(means);se=math.sqrt(variance(means)/C);half=critical*se
        rh=split_rhat(ps);change=max(abs(x-y)/max(x,1e-12) for x,y in zip(blockses,double));neffective=math.fsum(eff)
        quality=all(math.isfinite(x) for x in [value,se,half,rh,change,neffective]) and rh<=cfg['split_rhat_max'] and change<=cfg['blocking_relative_change_max'] and neffective>=cfg['effective_bins_min']
        covered=abs(value-ref)<=half
        for key,x,y in [('value',value,r['value']),('se',se,r['standard_error']),('halfwidth',half,r['half_width_95']),('rhat',rh,r['split_rhat']),('blocking',change,r['blocking_max_relative_change'])]:compare(key,x,y)
        assert quality==r['quality_pass'] and covered==r['covered'],r['id']
        rows.append({'id':r['id'],'value':value,'standard_error':se,'half_width_95':half,'covered':covered,'quality_pass':quality,'split_rhat':rh,'blocking_relative_change':change,'effective_bins':neffective})
    pilotseeds={s for xs in seeds['pilot'].values() for s in xs};assert allseeds.isdisjoint(pilotseeds)
    k=sum(r['covered'] for r in rows);pval=exact_pvalue(k,N,cfg['coverage_null']);cp=clopper_pearson(k,N,cfg['coverage_cp_confidence'])
    bias=avg([r['value']-ref for r in rows]);bias_se=math.sqrt(variance([r['value']-ref for r in rows])/N)
    bt=student_quantile((1.+cfg['bias_interval_confidence'])/2.,N-1);bias_ci=[bias-bt*bias_se,bias+bt*bias_se]
    assert k==recorded['covered_replicates'] and recorded['total_replicates']==N and recorded['raw_audit']=='PASS'
    for key,x,y in [('binomial_pvalue',pval,recorded['binomial_pvalue']),('mean_bias',bias,recorded['mean_bias']),('bias_se',bias_se,recorded['bias_standard_error'])]:compare(key,x,y)
    for j in range(2):
        compare('coverage_CP_endpoint',cp[j],recorded['coverage_clopper_pearson99'][j])
        compare('bias_interval_endpoint',bias_ci[j],recorded['mean_bias_interval99'][j])
    coverage=pval>=cfg['binomial_significance'] and cp[0]>=cfg['coverage_equivalence_lower'] and cp[1]<=cfg['coverage_equivalence_upper']
    bias_pass=bias_ci[0]>=-cfg['bias_equivalence_margin'] and bias_ci[1]<=cfg['bias_equivalence_margin']
    passed=coverage and bias_pass and all(r['quality_pass'] for r in rows)
    assert passed and recorded['status']=='PASS'
    report={'schema':'x06_independent_stdlib_raw_audit_v1','catalogue_id':'X06','status':'PASS',
        'raw_sha256':sha(a.raw),'recorded_summary_sha256':sha(a.recorded_summary),'auditor_sha256':sha(__file__),
        'implementation_independence':'stdlib only; pseudovalue means/SE, direct log-PMF binomial tails, inverted binomial CP limits, incomplete-beta Student CDF; no imports from frozen comparator',
        'replicates_audited':N,'chains_audited':N*C,'raw_paired_bins_audited':N*C*B,
        'all_seed_derivations_match':True,'pilot_final_seeds_disjoint':True,'all_denominator_bins_positive_16':all_positive_sign,
        'coverage_count':k,'observed_coverage':k/N,'binomial_two_sided_pvalue':pval,'coverage_CP99':cp,
        'mean_bias':bias,'bias_standard_error':bias_se,'bias_CI99':bias_ci,'student_t7':critical,'student_t399_995':bt,
        'quality_failures':[r['id'] for r in rows if not r['quality_pass']],
        'max_split_rhat':max(r['split_rhat'] for r in rows),'max_blocking_relative_change':max(r['blocking_relative_change'] for r in rows),
        'min_effective_bins':min(r['effective_bins'] for r in rows),'maximum_numeric_replay_differences':max_replay_error,
        'scope':'One declared finite-dtau interacting hopping fixture and actual eight-chain95% interval; no universal coverage claim or rare-tail calibration',
        'replicates':rows}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='replicates'},indent=2))
if __name__=='__main__':main()
