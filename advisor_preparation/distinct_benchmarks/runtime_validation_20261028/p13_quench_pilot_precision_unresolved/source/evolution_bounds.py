"""Prospective conditional P13 phase/log enclosures; never executed during preparation.

Uses unchanged P17 ball arithmetic; no numerical imports at module load.
No claim of directed-rounding interval arithmetic or a new numerical result.
"""
import math
from inherited_p17_certificate import RealBall,ComplexBall

def phase_exp_minus_i(a,values,time):
 """Enclose exp(-i*value*time), avoiding an assumed libm exp/cis accuracy."""
 np=a.np
 x=a.ball(np.asarray(values,dtype=a.ld))*a.ball(time)
 maximum=a.scalar_upper(np.max(np.abs(x.c)+x.r))
 exponent=0
 while maximum> a.ld('0.25'):
  maximum=maximum/2
  exponent+=1
  if exponent>32:raise RuntimeError('Predeclared phase scaling cap exceeded')
 y=ComplexBall(a.ball(np.zeros_like(x.c)),(-x)/a.ball(2**exponent))
 if np.max(np.abs(y.imag.c)+y.imag.r)>a.ld('0.25'):
  raise RuntimeError('Scaled phase magnitude not qualified')
 one=ComplexBall(a.ball(np.ones_like(x.c)),a.ball(np.zeros_like(x.c)))
 term=one;result=one
 for k in range(1,21):
  term=(term*y).divide_real(a.ball(k))
  result=result+term
 # e^(1/4)*(1/4)^21/21! < 2^-100. Each component is bounded by this tail.
 tail=(a.ball(1)/a.ball(2**100)).upper()
 result=ComplexBall(RealBall(a,result.real.c,a.inflate(result.real.r+tail)),RealBall(a,result.imag.c,a.inflate(result.imag.r+tail)))
 for _ in range(exponent):result=result*result
 return result

def _atanh_log_series(a,z):
 """2*atanh(z),96terms; |z|<=3/5 implies remainingtail<2^-120."""
 np=a.np
 if np.max(np.abs(z.c)+z.r)>a.ld('0.6'):
  raise RuntimeError('Log series interval outside registered3/5 radius')
 z2=z*z;power=z;value=a.ball(0)
 for k in range(96):
  value=value+power/a.ball(2*k+1)
  power=power*z2
 value=value*2
 tail=(a.ball(1)/a.ball(2**120)).upper()
 return RealBall(a,value.c,a.inflate(value.r+tail))

def log_positive(a,value):
 """Enclose ln(value) for a positive scalarball without relying on libm log."""
 if value.c.shape!=()or value.lower()<=0:
  raise RuntimeError('Strictlypositive scalarpurity bound required')
 # frexp selects an exactpower-of-two scaling only. A separate ballrange gate
 # validates that choice, so binary64 conversion is not an accuracy assumption.
 exponent=math.frexp(float(value.c))[1]
 if not -16<=exponent<=16:raise RuntimeError('Purity log scaling outside registeredrange')
 multiplier=2.0**(-exponent)
 scaled=value*a.ball(multiplier)
 if scaled.lower()<a.ld('0.25')or scaled.upper()>a.ld(2):
  raise RuntimeError('Log scaled interval outside[1/4,2]')
 z=(scaled-1)/(scaled+1)
 ln_m=_atanh_log_series(a,z)
 ln_2=_atanh_log_series(a,a.ball(1)/a.ball(3))
 return ln_m+ln_2*exponent

def spectral_action_ball(a,V,values,initial,time):
 """Actual recorded complete V/Lambda spectral formula, with phase/arithmeticball."""
 z=ComplexBall.exact(a,initial).reshape((len(values),1))
 # V is real for these real OBC Hubbard Hamiltonians; no finaleigenpair matching.
 if a.np.any(a.np.asarray(V).imag!=0):raise RuntimeError('Unexpected complex eigensystem')
 realV=a.ball(a.np.asarray(V).real)
 coeff=ComplexBall(realV.T.matmul(z.real),realV.T.matmul(z.imag))
 phase=phase_exp_minus_i(a,values,time).reshape((len(values),1))
 weighted=phase*coeff
 return ComplexBall(realV.matmul(weighted.real),realV.matmul(weighted.imag))

def eigensystem_operator_error(a,certificate,time):
 """Unitary functional calculus: alpha(s+1)+|time|eta for raw V vs physical H."""
 q=a.ball(certificate['alpha_upper'])*(a.ball(certificate['s_upper'])+1)
 return (q+a.ball(abs(time))*a.ball(certificate['eta_upper'])).upper()

def initial_ground_state_distance(a,certificate,values):
 """Conservative same-sector polar/Davis-Kahan bound; uniqueground only."""
 delta=a.ball(values[1])-a.ball(values[0])
 eta=a.ball(certificate['eta_upper']);gap=delta-eta*2
 if gap.lower()<=0:raise RuntimeError('Initial uniqueground gap not qualified')
 # Polarrecordedground column vs actualH ground. sqrt2 times projector bound.
 distance=(a.ball(2).sqrt()*(eta*2/gap)+a.ball(certificate['alpha_upper'])).upper()
 return {'gap_lower':gap.lower(),'recorded_ground_distance_upper':distance}

def entropy_from_purity_with_state_budget(a,purity,state_distance_upper):
 """|purity(rho)-purity(sigma)|<=4delta for normalized states; explicitlogball."""
 extra=(a.ball(state_distance_upper)*4).upper()
 ideal_purity=RealBall(a,purity.c,a.inflate(purity.r+extra))
 return -log_positive(a,ideal_purity)


def growth_rate_ball(a,entropy_balls):
 """Fixed nine-point OLS using actual rational coefficient balls."""
 if len(entropy_balls)!=33:raise RuntimeError('Complete33timeentropycurve required')
 result=a.ball(0)
 for k in range(4,13):
  weight=a.ball(2*(k-8))/a.ball(15)
  result=result+entropy_balls[k]*weight
 return result
