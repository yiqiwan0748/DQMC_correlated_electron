#!/usr/bin/env python3
"""Promote only the two fully audited weekend scientific questions; no job/solver."""
import argparse,copy,hashlib,json,tarfile
from datetime import datetime,timezone
from pathlib import Path

BASE='728f0a1642dfc4ac4c099478ec033f40981871c11101ee76bc7e299918d2b74c'
AUDIT14='285613cd21ac3c285eb8905cd2c83612da73726ed5f01606cabc915b6fd23300'
AUDIT19='63d6490dc1b9224ae3c1ef53be31abb56704d6fad242e55a32b7676753d5057b'
COMP14='fad0e1f41adb991b07a7d5cb9369d7ce2f1abf5e591804ce9f702ba92d7b2f17'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(v):return json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def load(p):return json.loads(Path(p).read_text())
def need(ok,msg):
    if not ok:raise ValueError(msg)
def write(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();db=a.repo_root.resolve()/'advisor_preparation/distinct_benchmarks'
    bp=db/'reporting/reports/a25_free_ring_flux_v16/ledger_100.json'
    need(sha(bp)==BASE,'Base ledger identity');base=load(bp)
    need(base['counts']==dict(PASS=47,FAIL=0,UNRESOLVED=2,NOT_RUN=51),'Base tally')
    packet=db/'weekend_external_20261025';mf=load(packet/'PACKET_MANIFEST.json')
    need(mf['status']=='ROOT_REVIEWED_SCIENTIFIC_PACKET','Packet not closed')
    data={};origins={};member_count=0
    for item in mf['auxiliary_files']:
        q=packet/item['path'];need(sha(q)==item['sha256'] and q.stat().st_size==item['bytes'],'Auxiliary identity')
    for arc in mf['archives']:
        q=packet/arc['path'];need(sha(q)==arc['sha256'] and q.stat().st_size==arc['bytes'],'Archive identity')
        with tarfile.open(q,'r:gz') as t:
            members=t.getmembers();need(len(members)==len(arc['members']),'Archive count')
            need(len({m.name for m in members})==len(members),'Duplicate member')
            for member,entry in zip(members,arc['members']):
                need(member.isfile() and not Path(member.name).is_absolute() and '..' not in Path(member.name).parts,'Unsafe archive member')
                b=t.extractfile(member).read();b.decode('utf-8')
                need(member.name==entry['public_member'] and hashlib.sha256(b).hexdigest()==entry['public_sha256'] and len(b)==entry['public_bytes'],'Member identity')
                key=entry['source_identity'];need(key not in data,'Duplicate scientific identity')
                data[key]=b;origins[key]=entry;member_count+=1
    def original(path,expected):
        need(path in data and origins[path]['original_sha256']==expected,'Original scientific binding: '+path)
        return json.loads(data[path])
    c14=original('results/x14_saved_compare_recovery_v1_pc1_20261024/comparison.json',COMP14)
    x14=original('results/x14_full_independent_saved_audit_v1_pc1_20261024/audit.json',AUDIT14)
    x19=original('results/x19_full_saved_audit_v4_pc1_20261024/audit.json',AUDIT19)
    c19=original('analysis_inputs/x19_full_packet_6466825_v1/results/comparison.json','5f5cd621630d4b2314ce8ba42ec230fd955c910f30ced991759d34b2842cf077')
    replay19=original('results/x19_full_saved_audit_v4_pc1_20261024/frozen_comparison_replay.json','5f5cd621630d4b2314ce8ba42ec230fd955c910f30ced991759d34b2842cf077')
    need(c19==replay19 and c19['status']==c19['catalogue_status']=='PASS' and c19['full_report_complete'] and not c19['pilot_excluded'],'Original full X19 comparison and replay')
    def gates(value):
        if isinstance(value,dict):
            found=[value] if {'passed','bound','relation'}.issubset(value) else []
            return found+[g for v in value.values() for g in gates(v)]
        if isinstance(value,list):return [g for v in value for g in gates(v)]
        return []
    g19=gates(c19)
    need(len(g19)==1864 and all(g['passed'] for g in g19),'All frozen X19 comparator gates')
    cancelled=original('results/x14_exact_sign_full_v3_20261023T2115Z/execution.json','7c687f96dcee9d5e52be95034f09716f53ce58f8823219c4c0b4a337a8a3134d')
    need(cancelled['scientific_completions']==0 and cancelled['status']=='UNRESOLVED','Original cancellation remains zero credit')
    need(c14['status']==c14['catalogue_status']=='PASS' and c14['mode']=='production' and not c14['pilot_excluded'],'Full X14 decision')
    need(len(c14['checks'])==10024 and all(v['passed'] for v in c14['checks']) and len(c14['details'])==140,'Full X14 gates')
    need(x14['status']==x14['scientific_status']==x14['independent_mathematics']=='PASS' and x14['grid_rows']==140,'Full X14 audit scope')
    need(x14['check_count']==len(x14['checks'])==260842 and all(v['passed'] for v in x14['checks']),'X14 audit decisions')
    need(x14['proof_statistics']['sign_requests']==33769 and x14['proof_statistics']['all_events_consumed'] and x14['full_census_replayed'] and x14['every_point_sign_independently_exact'],'X14 complete proof')
    need(len(x14['mandatory_full_witnesses'])==8 and all(v['eligible'] for v in x14['mandatory_full_witnesses']),'X14 physical shift witnesses')
    need(x19['status']==x19['full_mathematics']==x19['catalogue_status']=='PASS' and x19['registered_cases']==12 and x19['spin_sectors_per_case_per_backend']==5,'Full X19 audit scope')
    need(x19['check_count']==len(x19['checks'])==18462 and all(v['passed'] for v in x19['checks']),'X19 audit decisions')
    need(x19['frozen_comparator_gate_count']==1864 and x19['exact_replay_equal'] and x19['completed_saved_matrix_eigh_calls']==144,'X19 complete frozen replay')
    review=load(packet/'independent/promotion_eligibility_review.json')
    need(review['status']=='PASS_SCIENTIFIC_ELIGIBILITY_SUBJECT_TO_ROOT_PROVENANCE_PUBLIC_CLOSURE','Independent eligibility status')
    need(review['eligible_distinct_questions']==['X14','X19'] and review['maximum_credit_from_these_results']==2,'Independent eligibility')
    for name,item in review['local_reports'].items():need(sha(packet/'independent'/name)==item['sha256'],'Independent actual-result binding')
    old={r['id']:r for r in base['rows']};need(old['X14']['status']==old['X19']['status']=='NOT_RUN','New question identities')
    rows=copy.deepcopy(base['rows'])
    for row in rows:
        if row['id'] not in ['X14','X19']:continue
        row.update(status='PASS',reported_status='PASS',audit_issues=[],previous_public_row_sha256=hashlib.sha256(canon(old[row['id']])).hexdigest())
        if row['id']=='X14':
            row.update(reason='Full Lorentzian peak/census comparison of previously audited native and independent QuSpin Lehmann data passes all fixed gates, including shifted/overlapping visible-peak witnesses.',fixtures=140,comparisons=10024,
                metrics=dict(maximum_independent_full_curve_abs=x14['max_independent_full_curve_abs'],complete_exact_point_sign_requests=33769,distinct_exact_points=18267,mandatory_physical_witnesses=8),
                scope=dict(model='Previously audited four-site Hubbard U4 unique fixed(2,2) initial state; 8 channels ×3 branches ×5 broadenings plus20controls.',method='New spectrum-rendering/peak qualification adapters; no new eigenstate, eigensolver or sampler.',qualification='PADDED_BINARY64_OR_INDEPENDENT_EXACT_RATIONAL_POINT_SIGN_V2. Exact rational point signs; interval census and curvature/root enclosures conditional on disclosed binary64 arithmetic.',independent='Complete 140-row ISAAC saved-data audit plus local stdlib census/full-curve replay and direct Fraction verification of all18267unique points.',limit='One broadening-error question; no physical gap, algorithm speedup, formal directed-interval or whole-package claim.'),
                input_role='Completed native/reference renderings from cancelled6467601 plus new comparator6469838 and full independent audit6469951; original cancelled execution remains zero credit.',
                evidence=dict(packet_manifest='../../../weekend_external_20261025/PACKET_MANIFEST.json',original_full_audit_sha256=AUDIT14,full_comparison_sha256=COMP14,independent_compute_checks=260842),
                scientific_history=[dict(job='6466312',mode='Original failed reference pilot; retained public record',catalogue_credit=0),dict(job='6467268',mode='Exact-sign pilot; excluded',catalogue_credit=0),dict(job='6467601',mode='Cancelled original full process; completed renderings retained, no comparison',catalogue_credit=0),dict(job='6469838',mode='Separately admitted comparison of unchanged saved full data',catalogue_credit=1),dict(job='6469951',mode='Independent complete saved-data mathematical audit',catalogue_credit=0)])
        else:
            row.update(reason='Actual native finite-cluster geometry and generic ED outputs agree with independent QuSpin, full physical Fock reconstruction and analytical mixing/canonical controls for all twelve fixed cases.',fixtures=12,comparisons=1864,
                metrics=dict(local_independent_K_max_abs=review['x19']['max_local_K_abs'],local_independent_Fock_H_max_abs=review['x19']['max_local_H_abs'],complete_sector_blocks_both_backends=120,complete_eigenenergies_per_case_per_backend=1820),
                scope=dict(model='Matched open2x2, periodic2x2 and tilted4 two-orbital clusters; fixed totalN4, beta1, U0/4, t4=0/-.65; all5spinsectors.',method='Actual native geometry/generic ED compared with actual QuSpin; new observation and independent saved-matrix reconstruction adapters.',observables='Complete spectra, canonical spin-sector probabilities, physical spin/charge correlations, orbital coherence and transfer; own-geometry controls and cross-geometry/t4 contrasts.',runtime_repair='One explicitly prospectively pinned grp extension; all scientific functions/physical gates unchanged. Three failed audits preserved; historical failed process mappings are not recovered.',limit='One finite geometry/mixing question, without quantitative journal phase-diagram reproduction, interacting Lehmann spectrum, thermodynamic or whole-package claim.'),
                input_role='Original full6466825 and complete later saved-data audit6468775 with exact frozen comparator replay; local independent sparse-Fock/character check.',
                evidence=dict(packet_manifest='../../../weekend_external_20261025/PACKET_MANIFEST.json',original_full_audit_sha256=AUDIT19,independent_compute_checks=18462,independent_local_checks=4619),
                scientific_history=[dict(job='6466825',mode='Full fixed twelve-case scientific run',catalogue_credit=1),dict(job='6467037',mode='Failed saved-data audit v1; retained',catalogue_credit=0),dict(job='6467245',mode='Failed saved-data audit v2; retained',catalogue_credit=0),dict(job='6467556',mode='Failed saved-data audit v3; retained',catalogue_credit=0),dict(job='6468775',mode='Complete saved-data audit v4 and exact frozen comparator replay',catalogue_credit=0)])
    changed=['X14','X19'];preserved={}
    for row in rows:
        if row['id'] not in changed:
            need(canon(row)==canon(old[row['id']]),'Unaffected row changed')
            preserved[row['id']]=hashlib.sha256(canon(row)).hexdigest()
    ids={f'{p}{i:02d}' for p in ['A','J','P','X'] for i in range(1,26)}
    need(len(rows)==100 and {r['id'] for r in rows}==ids,'Exactly100distinct questions')
    counts={s:sum(r['status']==s for r in rows) for s in ['PASS','FAIL','UNRESOLVED','NOT_RUN']}
    need(counts==dict(PASS=49,FAIL=0,UNRESOLVED=2,NOT_RUN=49),'Counted repeated/unaudited idea')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    write(out/'ledger_100.json',dict(schema='public_weekend_reconciliation_leaf_v17',generated_utc=datetime.now(timezone.utc).isoformat(),tool_sha256=sha(__file__),counts=counts,category_counts={c:{s:sum(r['category']==c and r['status']==s for r in rows) for s in counts} for c in sorted({r['category'] for r in rows})},counting_rule=base['counting_rule'],scope_limit='Only X14 and X19 gain one question each. All98otherrows/all47priorPASS/bothUNRESOLVED unchanged; failures/cancellations/conditional arithmetic retained.',base_ledger='../a25_free_ring_flux_v16/ledger_100.json',base_ledger_sha256=BASE,base_preservation=dict(all98unaffected_rows_unchanged=True,all47prior_PASS_rows_unchanged=True,both_prior_UNRESOLVED_rows_unchanged=True,canonical_preserved_row_sha256=preserved),changed_ids=changed,packet_manifest='../../../weekend_external_20261025/PACKET_MANIFEST.json',packet_manifest_sha256=sha(packet/'PACKET_MANIFEST.json'),rows=rows))
    write(out/'audit.json',dict(schema='PUBLIC_WEEKEND_LEAF_CLOSURE_AUDIT_V17',status='PASS',ledger_sha256=sha(out/'ledger_100.json'),counts=counts,changed_ids=changed,source_data_archives=len(mf['archives']),regular_UTF8_members=member_count,all98other_rows_preserved=True,pilot_analysis_setup_credit=0,producer_or_new_scientific_job_executed=False))
    (out/'README.md').write_text('# Weekend results reconciled\n\n**49 PASS, 2 UNRESOLVED and 49 NOT_RUN** across exactly100 distinct catalogue questions. Only X14 and X19 change from the immutable [47-question checkpoint](../a25_free_ring_flux_v16/README.md). All98otherrow objects, all47previousPASS and both unresolved histories remain exact.\n\nX14 completes the fixed spectral-broadening/peak-shift comparison using unchanged saved native/QuSpin Lehmann data,10,024 gates and260,842 independent compute checks. The cancelled original process remains zero credit; its completed renderings are reused by a separately admitted recovery. Exact point signs and conditional binary64 interval-census qualifications remain distinct.\n\nX19 completes the twelve-case finite-geometry/orbital-mixing comparison against QuSpin and analytical controls,1,864 frozen gates and18,462 independent compute checks. Its three failed audits remain preserved beside the narrow, prospectively pinned runtime repair.\n\nThe [scientific evidence packet](../../../weekend_external_20261025/README.md) retains full saved arrays/proof streams, mathematical source, frozen gates, audit and failure histories. Later local saved-data checks are explicitly attributed and do not claim the original pinned runtime. Each idea counts once. No complete-package or thermodynamic reliability claim follows. [Closure audit](audit.json).\n'.replace('exactly100','exactly 100').replace('All98otherrow','All 98 other row').replace('all47previousPASS','all 47 previous PASS').replace(',10,',', 10,').replace('and260,','and 260,').replace(',1,864',', 1,864').replace('and18,','and 18,'))
    print(json.dumps(dict(status='PASS',counts=counts,ledger_sha256=sha(out/'ledger_100.json')),indent=2))

if __name__=='__main__':main()
