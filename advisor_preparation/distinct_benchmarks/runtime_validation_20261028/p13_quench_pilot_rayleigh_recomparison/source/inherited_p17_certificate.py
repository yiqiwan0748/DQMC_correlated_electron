"""NEW conditional extended-accumulation balls; unexecuted source candidate.

No numerical import or finfo call occurs at module import. All model/eigensolver
inputs remain represented binary64 values. Arithmetic is performed only by a
future separately admitted/pinned compute process, never preparation.
"""
from __future__ import annotations
import math

class Arithmetic:
    MAX_DIMENSION=4096
    MAX_DOT_LENGTH=4096
    def __init__(self,np):
        self.np=np
        self.ld=np.longdouble
        f=np.finfo(np.longdouble)
        if f.nmant<63:
            raise RuntimeError('Unqualified longdouble precision')
        if f.nmant<63 or f.eps!=self.ld(2)**(-int(f.nmant)):
            raise RuntimeError('Requires actual binary longdouble significand>=64bits')
        self.u=self.ld(f.eps)/self.ld(2)
        self.tiny=self.ld(f.tiny)
        self.calls={}
        self.metadata={'dtype':str(np.dtype(np.longdouble)),'nmant':int(f.nmant),
            'eps_decimal':str(f.eps),'unit_roundoff_decimal':str(self.u),
            'tiny_decimal':str(f.tiny),'arithmetic':'real longdouble; complex split into two real balls',
            'conditional_not_interval':True,'matmul_operation_allowance':'16k+32; positive-bound inflation128k+512',
            'maximum_dot_length':self.MAX_DOT_LENGTH}

    def check(self,x):
        a=self.np.asarray(x,dtype=self.ld)
        if not self.np.all(self.np.isfinite(a)):
            raise RuntimeError('Nonfinite arithmetic')
        if self.np.any((a!=0)&(self.np.abs(a)<self.tiny)):
            raise RuntimeError('Nonzero subnormal longdouble operand/result')
        return a

    def bump(self,name):
        self.calls[name]=self.calls.get(name,0)+1
        if sum(self.calls.values())>5000000:
            raise RuntimeError('Predeclared arithmetic call cap')

    def up(self,x):
        a=self.check(x)
        if self.np.any(a<0):
            raise RuntimeError('Upper envelope must be nonnegative')
        result=a.copy()
        nonzero=a!=0
        # Zero here is an already validated bound, a structural zero operand,
        # or a positive-envelope zero whose support is verified below.
        result[nonzero]=self.np.nextafter(a[nonzero],self.ld('inf'))
        return self.check(result)

    def down(self,x):
        a=self.check(x)
        result=a.copy()
        nonzero=a!=0
        # With mandatory underflow traps, an accepted exact floating zero is
        # structural or exact cancellation of represented operands. Moving it
        # to the smallest subnormal would itself violate the normality model.
        result[nonzero]=self.np.nextafter(a[nonzero],self.ld('-inf'))
        return self.check(result)

    def gamma(self,n):
        if not isinstance(n,int) or n<1 or n>1000000000:
            raise RuntimeError('Unregistered operation budget')
        nu=self.up(self.ld(n)*self.u)
        if nu>=self.ld('0.25'):
            raise RuntimeError('Invalid gamma regime')
        denominator=self.down(self.ld(1)-nu)
        return self.up(nu/denominator)

    def inflate(self,x,n=128):
        x=self.check(x)
        if self.np.any(x<0):
            raise RuntimeError('Nonnegative bound required')
        return self.up(x/self.down(self.ld(1)-self.gamma(n)))

    def positive_product(self,a,b):
        a=self.check(a);b=self.check(b)
        if self.np.any(a<0) or self.np.any(b<0):
            raise RuntimeError('Nonnegative envelope inputs')
        out=self.check(a*b)
        if self.np.any((out==0)&(a!=0)&(b!=0)):
            raise RuntimeError('Unproved zero in positive product')
        return out

    def positive_matmul(self,a,b):
        np=self.np;a=self.check(a);b=self.check(b)
        if a.ndim!=2 or b.ndim!=2 or a.shape[1]!=b.shape[0]:
            raise RuntimeError('Unregistered matrix shape')
        k=a.shape[1]
        if k<1 or k>self.MAX_DOT_LENGTH or max(a.shape+b.shape)>self.MAX_DIMENSION:
            raise RuntimeError('Fixed certificate dimension cap')
        if np.any(a<0) or np.any(b<0):
            raise RuntimeError('Positive product required')
        out=self.check(a@b)
        # uint16 exact support counts cannot overflow because k<=4096<65535.
        support=(a!=0).astype(np.uint16)@(b!=0).astype(np.uint16)
        if np.any((out==0)&(support!=0)) or np.any((out!=0)&(support==0)):
            raise RuntimeError('Unproved positive-matmul zero')
        return self.inflate(out,16*k+32)

    def norm_upper(self,center):
        np=self.np;a=np.abs(self.check(center)).ravel()
        if a.size==0:
            raise RuntimeError('Empty norm')
        scale=np.max(a)
        if scale==0:
            return self.ld(0) # every center entry was structurally zero
        q=self.check(a/scale)
        squares=self.positive_product(q,q)
        total=self.check(np.sum(squares,dtype=self.ld))
        if total==0:
            raise RuntimeError('Unproved positive norm zero')
        total=self.inflate(total,16*int(a.size)+64)
        result=self.check(scale*np.sqrt(total))
        return self.inflate(result,32)

    def scalar_upper(self,x):
        """Checked outward binary64 serialization for a NONNEGATIVE bound."""
        a=self.check(x)
        if a.shape!=() or a<0:
            raise RuntimeError('Expected nonnegative scalar upper bound')
        if a==0:
            return 0.0 # caller supplied a validated zero bound
        out=math.nextafter(float(a),math.inf)
        if not math.isfinite(out) or self.ld(out)<a:
            raise RuntimeError('Downward or nonfinite serialized upper bound')
        return out

    def ball(self,value,radius=0):
        return RealBall(self,value,radius)

class RealBall:
    def __init__(self,arithmetic,center,radius=0):
        self.a=arithmetic
        self.c=self.a.check(center)
        self.r=self.a.check(self.a.np.broadcast_to(radius,self.c.shape))
        if self.a.np.any(self.r<0):
            raise RuntimeError('Negative error radius')

    def _other(self,x):
        return x if isinstance(x,RealBall) else self.a.ball(x)

    def __add__(self,other):
        a=self.a;b=self._other(other);a.bump('add')
        c=a.check(self.c+b.c)
        radius=self.r+b.r+a.gamma(2)*(a.np.abs(self.c)+a.np.abs(b.c))
        return RealBall(a,c,a.inflate(radius,64))
    __radd__=__add__

    def __neg__(self):
        return RealBall(self.a,-self.c,self.r) # exact sign-bit change

    def __sub__(self,other):
        return self+-self._other(other)

    def __rsub__(self,other):
        return self._other(other)+-self

    def __mul__(self,other):
        a=self.a;b=self._other(other);a.bump('multiply')
        c=a.check(self.c*b.c)
        radius=a.positive_product(self.r,a.np.abs(b.c)+b.r)+a.positive_product(a.np.abs(self.c),b.r)
        radius+=a.gamma(2)*a.positive_product(a.np.abs(self.c),a.np.abs(b.c))
        return RealBall(a,c,a.inflate(radius,128))
    __rmul__=__mul__

    def __truediv__(self,other):
        a=self.a;b=self._other(other);a.bump('divide')
        denominator=a.down(a.np.abs(b.c)-b.r)
        if a.np.any(denominator<=0):
            raise RuntimeError('Division interval contains zero')
        c=a.check(self.c/b.c)
        radius=(self.r+a.positive_product(a.np.abs(c),b.r)+
            a.gamma(2)*a.positive_product(a.np.abs(c),a.np.abs(b.c)))/denominator
        return RealBall(a,c,a.inflate(radius,128))

    def matmul(self,other):
        a=self.a;b=self._other(other);a.bump('matmul')
        if self.c.ndim!=2 or b.c.ndim!=2 or self.c.shape[1]!=b.c.shape[0]:
            raise RuntimeError('Expected explicit2D matrix product')
        k=self.c.shape[1]
        c=a.check(self.c@b.c)
        magnitude=a.positive_matmul(a.np.abs(self.c),a.np.abs(b.c))
        radius=a.positive_matmul(a.np.abs(self.c),b.r)+a.positive_matmul(self.r,a.np.abs(b.c))+a.positive_matmul(self.r,b.r)
        radius+=a.gamma(16*k+32)*magnitude
        return RealBall(a,c,a.inflate(radius,128*k+512))

    def sum(self,axis=None):
        a=self.a;a.bump('sum');np=a.np
        n=int(self.c.size if axis is None else self.c.shape[axis])
        c=a.check(np.sum(self.c,axis=axis,dtype=a.ld))
        magnitude=np.sum(np.abs(self.c),axis=axis,dtype=a.ld)
        radius=np.sum(self.r,axis=axis,dtype=a.ld)+a.gamma(16*n+32)*a.inflate(magnitude,16*n+32)
        return RealBall(a,c,a.inflate(radius,64*n+128))

    def sqrt(self):
        a=self.a;lower=a.down(self.c-self.r)
        if a.np.any(lower<=0):
            raise RuntimeError('Sqrt requires strictly positive lower bound')
        c=a.check(a.np.sqrt(self.c))
        denominator=a.down(a.np.sqrt(lower)+c)
        radius=self.r/denominator+a.gamma(4)*a.np.abs(c)
        return RealBall(a,c,a.inflate(radius,128))

    def norm_upper(self):
        return self.a.inflate(self.a.norm_upper(self.c)+self.a.norm_upper(self.r),64)

    def upper(self):
        return self.a.up(self.c+self.r)

    def lower(self):
        return self.a.down(self.c-self.r)

    def reshape(self,shape):
        return RealBall(self.a,self.c.reshape(shape),self.r.reshape(shape))

    @property
    def T(self):
        return RealBall(self.a,self.c.T,self.r.T)

    def __getitem__(self,key):
        return RealBall(self.a,self.c[key],self.r[key])

class ComplexBall:
    def __init__(self,real,imag):
        if real.c.shape!=imag.c.shape or real.a is not imag.a:
            raise RuntimeError('Split-complex shape/arithmetic mismatch')
        self.real=real;self.imag=imag;self.a=real.a

    @classmethod
    def exact(cls,a,z):
        z=a.np.asarray(z)
        return cls(a.ball(z.real),a.ball(z.imag))

    def __add__(self,other):
        return ComplexBall(self.real+other.real,self.imag+other.imag)

    def __sub__(self,other):
        return ComplexBall(self.real-other.real,self.imag-other.imag)

    def __mul__(self,other):
        if isinstance(other,RealBall) or not isinstance(other,ComplexBall):
            return ComplexBall(self.real*other,self.imag*other)
        return ComplexBall(self.real*other.real-self.imag*other.imag,
                           self.real*other.imag+self.imag*other.real)
    __rmul__=__mul__

    def divide_real(self,denominator):
        return ComplexBall(self.real/denominator,self.imag/denominator)

    def matmul(self,other):
        return ComplexBall(self.real.matmul(other.real)-self.imag.matmul(other.imag),
                           self.real.matmul(other.imag)+self.imag.matmul(other.real))

    def conjugate(self):
        return ComplexBall(self.real,-self.imag)

    def sum(self,axis=None):
        return ComplexBall(self.real.sum(axis),self.imag.sum(axis))

    def norm_squared(self):
        return (self.real*self.real+self.imag*self.imag).sum()

    def norm_upper(self):
        # Triangle upper bound avoids an unaccounted complex absolute/sqrt.
        return self.a.inflate(self.real.norm_upper()+self.imag.norm_upper(),64)

    def reshape(self,shape):
        return ComplexBall(self.real.reshape(shape),self.imag.reshape(shape))

    def __getitem__(self,key):
        return ComplexBall(self.real[key],self.imag[key])


def eigensystem_certificate(a,H,V,values):
    """H is a real ball for the exact physical model; V/Lambda are actual binary64."""
    np=a.np
    if H.c.shape!=V.shape or V.shape[0]!=V.shape[1] or len(values)!=V.shape[0]:
        raise RuntimeError('Complete square eigensystem required')
    n=len(values)
    VV=a.ball(V)
    lam=a.ball(np.asarray(values,dtype=a.ld).reshape((1,n)))
    gram=VV.T.matmul(VV)-a.ball(np.eye(n,dtype=a.ld))
    e=gram.norm_upper()
    if not 0<=e<1:
        raise RuntimeError('Polar Gram error must be<1')
    beta=(H.matmul(VV)-VV*lam).norm_upper()
    one=a.ball(1)
    eb=a.ball(e)
    s=(one+eb).sqrt()
    alpha=eb/(one+(one-eb).sqrt())
    normH=H.norm_upper()
    normLambda=a.up(np.max(np.abs(np.asarray(values,dtype=a.ld))))
    eta=(a.ball(normH)*eb+a.ball(beta)*s+a.ball(normLambda)*alpha*(one+s)).upper()
    return {'dimension':n,'gram_error_upper':e,'residual_frobenius_upper':beta,
            'eta_upper':eta,'alpha_upper':alpha.upper(),'s_upper':s.upper(),
            'H_norm_upper':normH,'Lambda_norm_upper':normLambda}


def spectral_response(a,H,V,values,diagA,diagO,certificate):
    """Complete recorded-V spectral sum plus H/polar/arithmetic radii."""
    np=a.np;n=len(values)
    if n<2:
        raise RuntimeError('Response requires complete excited-state complement')
    VV=a.ball(V)
    v0=VV[:,0].reshape((1,n))
    excited=VV[:,1:]
    A=a.ball(np.asarray(diagA,dtype=a.ld).reshape((1,n)))
    O=a.ball(np.asarray(diagO,dtype=a.ld).reshape((1,n)))
    left=(v0*A).matmul(excited)
    right=(v0*O).matmul(excited)
    denom=a.ball(np.asarray(values[1:],dtype=a.ld))-a.ball(values[0])
    response=((left*right)/denom).sum()*2
    Delta=a.ball(values[1])-a.ball(values[0])
    eta=a.ball(certificate['eta_upper'])
    g=Delta-eta*2
    if g.lower()<=0:
        raise RuntimeError('No qualified same-sector response gap')
    p=eta*2/g
    # Inverse minimum denominators are rounded down before forming reciprocal balls.
    c=a.ball(1)/a.ball(min(a.ld(1),g.lower()))
    cp=a.ball(1)/a.ball(min(a.ld(1),Delta.lower()))
    BH=(p*2/g+c*cp*(eta*2+p)+p)*a.ld('0.5') #2||A||||O||=.5
    q=a.ball(certificate['alpha_upper'])*(a.ball(certificate['s_upper'])+1)
    BV=(q*(a.ball(certificate['gram_error_upper'])+2)/Delta)*a.ld('0.5')
    serialized_center=float(response.c)
    # The stored center is binary64. Its conversion error is propagated when
    # this receipt is used as the h0 response in the field-window comparison.
    serialized_difference=a.ball(serialized_center)-response
    serialization_and_arithmetic=a.inflate(a.np.abs(serialized_difference.c)+serialized_difference.r,128)
    total=(a.ball(BH.upper())+a.ball(BV.upper())+a.ball(serialization_and_arithmetic)).upper()
    return {'center':serialized_center,'radius':a.scalar_upper(total),
            'arithmetic_radius':a.scalar_upper(serialization_and_arithmetic),
            'H_perturbation_radius':a.scalar_upper(BH.upper()),
            'recorded_V_radius':a.scalar_upper(BV.upper())}
