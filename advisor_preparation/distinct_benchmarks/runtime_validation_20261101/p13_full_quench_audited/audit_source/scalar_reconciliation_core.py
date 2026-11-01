"""SOURCE ONLY: independent saved-artifact scalar reconciliation, not a solver.

No main, filesystem, package imports, host operations or automatic credit.
The separately reviewed PC entry must qualify this source and its own complete
stdlib/loaded closure BEFORE invoking a data loader or these numeric functions.
"""
from fractions import Fraction
import math

TIMES=tuple(range(33));BACKENDS=('native','tenpy');ARMS=('coarse','primary')
LIMITS={'matrix':1e-11,'gap':1e-3,'norm':1e-12,'initial_residual':1e-10,'spectral_identity':1e-10,
        'conversion':1e-10,'library_entropy':1e-10,'primary_radius':4e-8,'cross_entropy':1e-7,
        'cross_purity':1e-8,'slope':1e-6,'arm_entropy':1e-7,'curve_range':1e-4}
COVERAGE={'model_certificates':4,'initial_state_certificates':2,'state_certificates':132,
 'times_per_backend_arm':33,'primary_cross_backend_S2_checks':33,'primary_cross_backend_purity_checks':33,
 'primary_coarse_curve_checks':66,'growth_rate_certificates':4,'growth_window_points':9,'primary_growth_slope_checks':1}

class AuditHold(ValueError):pass

def need(ok,message):
 if not ok:raise AuditHold(message)

def rational(value):
 need(type(value)in(int,float),'Finite actual int/binary64 required')
 if type(value)is float:need(math.isfinite(value),'Nonfinite binary64 artifact');return Fraction.from_float(value)
 return Fraction(value)

def decimal_rational(value):
 need(type(value)is str and 0<len(value)<=128 and all(c in'0123456789+-.eE'for c in value),'Bounded decimal scalar required')
 parts=value.lower().split('e');need(len(parts)<=2,'Single decimal exponent required')
 if len(parts)==2:
  exponent=parts[1];digits=exponent.lstrip('+-');need(bool(digits)and digits.isascii()and digits.isdigit()and len(digits)<=6,'Bounded decimal exponent syntax required');need(abs(int(exponent))<=5000,'Decimal exponent outside fixed5000 cap')
 try:return Fraction(value)
 except (ValueError,ZeroDivisionError,OverflowError)as error:raise AuditHold('Invalid decimal scalar')from error

def ball(value):
 need(type(value)is dict and set(value)=={'center','radius'},'Exact serialized ball required')
 c=rational(value['center']);r=rational(value['radius']);need(r>=0,'Negative radius');return c,r

def envelope_difference(left,right):
 c,r=ball(left);d,s=ball(right);return abs(c-d)+r+s

def exact_ols_point(curve):
 """Exact rational OLS of actual binary64 library outputs, intercept included.

Does not call or reproduce the original longdouble ball implementation. Its
serialized output radius is not used as a substituted input error budget.
 """
 need(set(curve)==set(TIMES),'All33 frozen observation times required')
 return sum((Fraction(2*(k-8),15)*rational(curve[k])for k in range(4,13)),Fraction(0))

def check_runtime_declared(value,role,stage):
 need(type(value)is dict and value.get('status')=='PASS','Original runtime guard not PASS')
 need(value.get('scientific_completions')==0 if role=='native'else value.get('scientific_credit')==0,'Runtime records cannot carry scientific credit/completions')
 # Source-native Julia and Python guards have different schemas and nested maps.
 # The complete original guard is preserved/hash-bound; a separate metadata
 # auditor must validate its closure and same-PID/raw-map/file associations.
 return {'role':role,'stage':stage,'status':'DECLARED_PASS_REQUIRES_ORIGINAL_CLOSURE_JOIN'}

def check_producer(raw,backend,report_by_key):
 schemas={'native':'P13_NATIVE_MPO_ENTROPY_DENSE_ADAPTER_CORE_V2','tenpy':'P13_TENPY_EXACTDIAG_ENTROPY_PRODUCER_CORE_V2'}
 need(raw.get('schema')==schemas[backend]and raw.get('status')=='PRODUCER_COMPLETE_NOT_A_SCIENTIFIC_PASS','Exact completed producer schema required')
 need(type(raw.get('scientific_credit'))is int and raw['scientific_credit']==0 and raw.get('pilot_only')is False,'Full producer data only, zero credit')
 coords=[list((u,d))for u in range(64)if u.bit_count()==3 for d in range(64)if d.bit_count()==3]
 # These coordinates are exact integer metadata, not an eigensystem/physics test.
 need(raw.get('coordinates')==coords,'Exact balanced400 canonical coordinate order differs')
 need(type(raw.get('models'))is list and len(raw['models'])==2 and [m.get('U')for m in raw['models']]==[8.0,4.0],'Both original Ui8/Uf4 models required')
 for model in raw['models']:
  for name in('Hamiltonian','eigenvectors'):
   pack=model[name];need(set(pack)=={'shape','order','real','imag'}and pack['shape']==[400,400]and pack['order']in('F','C')and len(pack['real'])==len(pack['imag'])==160000,'Complete400x400 raw model pack required')
   need(all(type(v)in(int,float)and math.isfinite(v)for v in pack['real'])and all(type(v)in(int,float)and math.isfinite(v)and v==0 for v in pack['imag']),'Finite real original eigensystem/matrix required')
  values=model['eigenvalues'];need(len(values)==400 and all(type(v)in(int,float)and math.isfinite(v)for v in values)and all(a<=b for a,b in zip(values,values[1:])),'Complete finite sorted400 eigenvalues required')
 rows=raw.get('samples');need(type(rows)is list and len(rows)==66,'All66 original samples required');seen={}
 for row in rows:
  k=row['time_numerator'];arm=row['arm'];need(type(k)is int and k in TIMES and arm in ARMS and type(row['time_denominator'])is int and row['time_denominator']==8 and rational(row['time'])==Fraction(k,8),'Exact frozen time/arm differs')
  key=(backend,arm,k);need(key not in seen,'Repeated producer sample key');seen[key]=row
  cert=report_by_key[key];need(rational(row['library_S2'])==rational(cert['library_S2_original'])and rational(row['raw_norm'])==rational(cert['raw_norm_reported_original']),'Original library entropy/raw norm not associated with certificate')
  if backend=='native':
   controls={'coarse':(128,1e-24),'primary':(256,1e-28)};budget,cut=controls[arm]
   need(row['maxdim']==budget and row['cutoff']==cut and row['actual_normalize_call']=='LinearAlgebra.normalize!(psi) before actual user entropy','Native registered controls/actual call differ')
   need(row['actual_normalize_call']==cert['actual_normalize_call']and row['pre_entropy_MPS_norm_before_normalize']==cert['pre_entropy_MPS_norm_before_normalize_original']and row['pre_entropy_MPS_norm_after_normalize']==cert['pre_entropy_MPS_norm_after_normalize_original'],'Native original normalization diagnostics differ')
  else:
   controls={'coarse':(128,1e-12),'primary':(256,1e-14)};budget,cut=controls[arm]
   need(row['external_bond_budget']==budget and row['from_full_cutoff']==cut,'TeNPy actual registered controls differ')
  state=row['state'];need(state['axes']==['left','right','physical']and len(state['sites'])==6 and state['export_stage']==('AFTER_ACTUAL_USER_ENTROPY_CALL'if backend=='native'else'AFTER_ACTUAL_TENPY_ENTROPY_CALL'),'Complete actual post-entropy MPS stage required')
  amplitudes=state['complete_site_amplitudes'];need(set(amplitudes)=={'shape','order','real','imag'}and amplitudes['shape']==[4096]and amplitudes['order']in('F','C')and len(amplitudes['real'])==len(amplitudes['imag'])==4096,'Original complete4096 post-MPS export required')
  need(all(type(v)in(int,float)and math.isfinite(v)for v in amplitudes['real'])and all(type(v)in(int,float)and math.isfinite(v)for v in amplitudes['imag']),'Nonfinite complete post-MPS export')
  bonds=[];previous=1
  for site in state['sites']:
   shape=site['shape'];need(len(shape)==3 and all(type(v)is int and v>0 for v in shape)and shape[0]==previous and shape[2]==4 and shape[0]<=64 and shape[1]<=64,'MPS original tensor shape/bond continuity differs');previous=shape[1];bonds.append(shape[1])
   size=shape[0]*shape[1]*shape[2];need(set(site)=={'shape','order','real','imag'}and site['order']in('F','C')and len(site['real'])==len(site['imag'])==size,'Complete original MPS tensor pack required')
   need(all(type(v)in(int,float)and math.isfinite(v)for v in site['real'])and all(type(v)in(int,float)and math.isfinite(v)for v in site['imag']),'Nonfinite MPS tensor')
  need(previous==1 and bonds[:-1]==cert['bond_ranks_original']and row['observed_max_bond']==max(bonds[:-1])and max(bonds[:-1])<=64 and max(bonds[:-1])<budget,'Actual rank/certificate/budget differs')
  for name in('actual_raw_evolved_state','actual_normalized_MPS_input_sector_state'):
   p=row[name];need(set(p)=={'shape','order','real','imag'}and p['shape']==[400]and p['order']in('F','C')and len(p['real'])==len(p['imag'])==400 and all(type(v)in(int,float)and math.isfinite(v)for v in p['real'])and all(type(v)in(int,float)and math.isfinite(v)for v in p['imag']),'Complete original400 state pack required')
 need(set(seen)=={(backend,a,k)for a in ARMS for k in TIMES},'Producer exact33x2 coverage differs')
 for k in TIMES:
  left=seen[backend,'coarse',k];right=seen[backend,'primary',k]
  need(all(left[n]==right[n]for n in('actual_raw_evolved_state','raw_norm','actual_normalized_MPS_input_sector_state')),'Arms no longer share actual original spectral state')
 for stage in('before','after'):check_runtime_declared(raw['runtime_'+stage],backend,stage)
 return {'backend':backend,'models':2,'samples':66,'actual_library_and_normalization_associations':66,'checked_pack_scope':['Hamiltonian','eigenvectors','eigenvalues','sample_raw_and_normalized400states','sample_complete4096amplitudes','sample_MPS_site_tensors'],'not_checked_pack_scope':['actual_initial_state','represented_MPO','local_operator_tables','TeNPy_pipe_phase_exports','Schmidt_value_exports'],'numerical_eigensystem_or_tensor_recompute':False}

def reconcile_report(report,expected_gate_schedule):
 """Requires the caller to have passed new own-PC source/runtime admission.

This is an independent exact-scalar audit of saved certificates, not a second
implementation of full-MPO/eigensystem/full-tensor/longdouble certification.
 """
 need(report.get('schema')=='P13_FULL_STATE_ENTROPY_COMPARISON_CORE_V1'and report.get('status')=='PASS_FULL_PENDING_INDEPENDENT_AUDIT','Only exact completed full comparison may enter this success audit')
 need(report.get('catalogue_id')=='P13'and report.get('catalogue_status')=='NOT_RUN'and type(report.get('scientific_credit'))is int and report['scientific_credit']==0 and report.get('pilot_only')is False and report.get('comparison_run_completed')is True,'Original full result/zero-credit invariants differ')
 need(report.get('conditional_floating_certificate_not_interval_arithmetic')is True and report.get('no_automatic_catalogue_credit')is True,'Original conditional/no-credit scope differs')
 need(report.get('failed_gate_names')==[]and not any(k in report for k in('error','after_guard_error')),'Original errors/failures cannot be erased')
 need(report.get('coverage')==report.get('expected_coverage')==COVERAGE,'Exact full frozen comparison coverage differs')
 gates=report['gates'];need(type(gates)is list and len(gates)==len(expected_gate_schedule)==11692,'All ordered frozen gate occurrences required')
 for actual,expected in zip(gates,expected_gate_schedule):
  need(actual['name']==expected['name']and actual['kind']==expected['kind']and actual['status']=='PASS','Ordered frozen gate occurrence/status differs')
  for field in('limit','required_minimum','required_lower'):
   need((field in actual)==(field in expected),'Frozen gate threshold field differs')
   if field in expected:need(rational(actual[field])==rational(expected[field]),'Frozen threshold changed')
  if 'limit'in expected:
   value=rational(actual['value_upper']);need(0<=value<=rational(expected['limit']),'Stronger serialized conservative gate envelope not within frozen threshold; audit HOLD, no radius shrinking')
   need(0<=decimal_rational(actual['qualified_longdouble_value'])<=rational(expected['limit']),'Original recorded longdouble predicate differs')
  if 'required_minimum'in expected:need(decimal_rational(actual['gap_lower_decimal'])>=rational(expected['required_minimum']),'Original fixed same-sector gap not qualified')
  if 'required_lower'in expected:need(decimal_rational(actual['range_lower_decimal'])>=rational(expected['required_lower']),'Original fixed full-curve range not exercised')
 need(len(report['model_certificates'])==4 and {(m['backend'],m['U'])for m in report['model_certificates']}=={(b,u)for b in BACKENDS for u in(8.0,4.0)},'Four complete model certificate keys required')
 for m in report['model_certificates']:
  need(m['dimension']==m['certificate']['dimension']==400,'Model certificate dimension differs')
  for n in('matrix_error_upper','actual_MPO_error_upper'):need(0<=rational(m[n])<=rational(LIMITS['matrix']),'Frozen matrix/MPO error threshold differs')
  for n,v in m['certificate'].items():
   if n!='dimension':need(rational(v)>=0,'Negative original model certificate upper bound')
 need(len(report['initial_state_certificates'])==2 and {r['backend']for r in report['initial_state_certificates']}==set(BACKENDS),'Both same-sector initial-ray certificates required')
 initial={r['backend']:r for r in report['initial_state_certificates']}
 for row in initial.values():
  need(0<=rational(row['initial_residual_upper'])<=rational(LIMITS['initial_residual']),'Initial residual upper gate differs')
  need(decimal_rational(row['Rayleigh_excited_separation_lower_decimal'])>0,'Current repaired Rayleigh separation nonpositive')
  ic,ir=ball(row['Rayleigh_energy_imaginary_ball']);need(ic-ir<=0<=ic+ir,'Hermitian Rayleigh imaginary enclosure excludes zero')
  need(rational(row['initial_physical_state_distance_upper'])>=0 and rational(row['actual_initial_norm_upper'])>0,'Actual repaired ray/normalization bound invalid')
 states=report['state_certificates'];need(type(states)is list and len(states)==132,'All132 full state certificates required');by={}
 for row in states:
  key=(row['backend'],row['arm'],row['time_numerator']);need(key not in by and row['backend']in BACKENDS and row['arm']in ARMS and type(row['time_numerator'])is int and row['time_numerator']in TIMES and type(row['time_denominator'])is int and row['time_denominator']==8,'Unique exact state certificate key required');by[key]=row
  for n in('raw_evolved_norm_squared','normalized_input_norm_squared','post_MPS_raw_norm_squared','Schmidt_probability_sum'):
   c,r=ball(row[n]);need(c-r>0 and abs(c-1)+r<=rational(LIMITS['norm']),'Serialized full norm/probability envelope exceeds frozen norm gate')
  need(row['actual_initial_norm_upper_multiplied_into_operator_error']==initial[row['backend']]['actual_initial_norm_upper']and row['initial_ground_and_column_normalization_distance_upper']==initial[row['backend']]['initial_physical_state_distance_upper'],'Current repaired initial-ray/actual-norm association differs')
  for field,limit in [('actual_spectral_state_identity_error_upper','spectral_identity'),('normalization_formula_identity_error_upper','conversion'),('complete_MPS_conversion_error_upper','conversion'),('complete_tensor_export_identity_error_upper','conversion'),('library_own_state_entropy_discrepancy_plus_budgets_upper','library_entropy')]:need(0<=rational(row[field])<=rational(LIMITS[limit]),'Serialized state error threshold differs: '+field)
  need(rational(row['outside_balanced_sector_weight_upper_retained'])>=0,'Outside-sector leakage erased/negative')
  for field in('fresh_normalized_MPS_purity','physical_purity'):
   c,r=ball(row[field]);need(c-r>0,'Nonpositive purity lower endpoint')
  c,r=ball(row['physical_entropy_centered_on_actual_library_value']);need(c==rational(row['library_S2_original']),'Physical entropy not centered on original actual library value')
  if row['arm']=='primary':need(r<=rational(LIMITS['primary_radius']),'Frozen primary entropy radius exceeded')
  for allowance in('TeNPy_omitted_purity_allowance','TeNPy_omitted_entropy_allowance','TeNPy_Schmidt_normalization_entropy_allowance'):need(rational(row[allowance])>=0,'Negative omission/normalization allowance')
  ownc,ownr=ball(row['fresh_normalized_MPS_entropy']);need(abs(c-ownc)+ownr+rational(row['TeNPy_omitted_entropy_allowance'])+rational(row['TeNPy_Schmidt_normalization_entropy_allowance'])<=rational(LIMITS['library_entropy']),'Stronger serialized library/own-state entropy gate exceeds threshold')
  if row['backend']=='native':need(row['TeNPy_omitted_purity_allowance']==row['TeNPy_omitted_entropy_allowance']==row['TeNPy_Schmidt_normalization_entropy_allowance']==0,'Native record acquired TeNPy omission budget')
 need(set(by)=={(b,a,k)for b in BACKENDS for a in ARMS for k in TIMES},'Exact132 full state key coverage differs')
 scalar_checks={'cross_backend_S2':0,'cross_backend_purity':0,'primary_coarse_S2':0,'growth_point_enclosures':0,'primary_slope':0,'curve_range':0}
 for k in TIMES:
  left=by['native','primary',k];right=by['tenpy','primary',k]
  need(envelope_difference(left['physical_entropy_centered_on_actual_library_value'],right['physical_entropy_centered_on_actual_library_value'])<=rational(LIMITS['cross_entropy']),'Stronger serialized cross-backend S2 bound fails');scalar_checks['cross_backend_S2']+=1
  need(envelope_difference(left['physical_purity'],right['physical_purity'])<=rational(LIMITS['cross_purity']),'Stronger serialized cross-backend purity bound fails');scalar_checks['cross_backend_purity']+=1
  for backend in BACKENDS:
   need(envelope_difference(by[backend,'primary',k]['physical_entropy_centered_on_actual_library_value'],by[backend,'coarse',k]['physical_entropy_centered_on_actual_library_value'])<=rational(LIMITS['arm_entropy']),'Stronger serialized arm-stability bound fails');scalar_checks['primary_coarse_S2']+=1
 fits=report['growth_rates'];need(type(fits)is list and len(fits)==4,'Allfour fixed nine-point fits required');growth={};fit_evidence=[]
 for row in fits:
  key=(row['backend'],row['arm']);need(key in{(b,a)for b in BACKENDS for a in ARMS}and key not in growth and row['time_numerators']==list(range(4,13))and row['weight_numerators']==[-8,-6,-4,-2,0,2,4,6,8]and type(row['weight_denominator'])is int and row['weight_denominator']==15,'Original OLS window/weights/key differs')
  growth[key]=row;point=exact_ols_point({k:by[key[0],key[1],k]['library_S2_original']for k in TIMES});center,radius=ball({'center':row['center'],'radius':row['radius']});difference=abs(point-center);need(difference<=radius,'Exact original-output OLS point outside saved fit enclosure');scalar_checks['growth_point_enclosures']+=1
  fit_evidence.append({'backend':key[0],'arm':key[1],'exact_rational_point':str(point),'serialized_center_error_exact':str(difference),'saved_radius_exact':str(radius),'scope':'Exact point OLS of actual binary64 library outputs; internal longdouble radius construction is source-bound and not independently reproduced.'})
 need(set(growth)=={(b,a)for b in BACKENDS for a in ARMS}and report['growth_rate_status']=='COMPLETE_FIXED_NINE_POINT_OLS','Four exact fit keys/status differ')
 left=growth['native','primary'];right=growth['tenpy','primary'];need(envelope_difference({'center':left['center'],'radius':left['radius']},{'center':right['center'],'radius':right['radius']})<=rational(LIMITS['slope']),'Stronger serialized primary growth-slope bound fails');scalar_checks['primary_slope']=1
 ranges=[]
 for backend in BACKENDS:
  balls=[ball(by[backend,'primary',k]['physical_entropy_centered_on_actual_library_value'])for k in TIMES];lower=max(c-r for c,r in balls)-min(c+r for c,r in balls);need(lower>=rational(LIMITS['curve_range']),'Stronger serialized complete-curve range not exercised');scalar_checks['curve_range']+=1;ranges.append({'backend':backend,'exact_serialized_envelope_range_lower':str(lower)})
 need(report['nontrivial_full_curve_status']=='CHECKED_FIXED_FULL_CURVE_NO_REPLACEMENT','Full curve status differs')
 for stage in('before','after'):check_runtime_declared(report['runtime_'+stage],'comparison',stage)
 return {'status':'PASS_INDEPENDENT_EXACT_SCALAR_RECONCILIATION_ONLY_PENDING_WHOLE_AUDIT','scientific_credit':0,'catalogue_status':'NOT_RUN','ordered_gate_occurrences':11692,'states':132,'scalar_checks':scalar_checks,'growth_fits':fit_evidence,'full_ranges':ranges,'original_conditional_longdouble_certificates_recomputed':False,'full_H_MPO_eigensystem_tensor_purity_log_certificates_recomputed':False,'source_bound_certificates_and_runtime_still_need_independent_metadata_and_whole_artifact_join':True},by
