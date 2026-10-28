"""P13 source-only exact combinatorics and fixed growth weights; no numerical imports."""
from fractions import Fraction
from itertools import product
L=6
NUP=NDN=3
CUT=3
U_INITIAL=8.0
U_FINAL=4.0
TIMES=tuple(Fraction(k,8)for k in range(33))
PILOT_TIMES=TIMES[:5]
GROWTH_INDICES=tuple(range(4,13))
GROWTH_WEIGHTS=tuple(Fraction(2*(k-8),15)for k in GROWTH_INDICES)
ARMS={'coarse':{'native_maxdim':128,'native_cutoff':1e-24,'tenpy_budget':128,'tenpy_cutoff':1e-12},'primary':{'native_maxdim':256,'native_cutoff':1e-28,'tenpy_budget':256,'tenpy_cutoff':1e-14}}
def coordinates():
 return tuple((u,d)for u in range(1<<L)if u.bit_count()==NUP for d in range(1<<L)if d.bit_count()==NDN)
def regroup_sign(u,d):
 return -1 if sum(((d>>j)&1)*sum((u>>k)&1 for k in range(j+1,L))for j in range(L))%2 else 1
def canonical_site_digits(u,d):
 return tuple(((u>>j)&1)+2*((d>>j)&1)for j in range(L))
def canonical_site_linear_index(u,d):
 return sum(x*4**j for j,x in enumerate(canonical_site_digits(u,d)))
def exact_growth_rate(values):
 if len(values)!=33:raise ValueError('Complete33-time curve required')
 if any(type(x)not in(int,Fraction)for x in values):raise ValueError('This exact helper accepts rational metadata controls only; production uses growth_rate_ball')
 return sum(w*values[k]for w,k in zip(GROWTH_WEIGHTS,GROWTH_INDICES))
def growth_absolute_error_bound(radii):
 if len(radii)!=33:raise ValueError('Complete33-time uncertainty curve required')
 return sum(abs(w)*radii[k]for w,k in zip(GROWTH_WEIGHTS,GROWTH_INDICES))
