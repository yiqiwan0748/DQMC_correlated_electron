# Independent analysis of the retained failed X14 root bracket

This is a new, post-failure diagnostic of saved trace data. It cannot qualify
X14, repair its reference renderer, change a gate or grant scientific credit.
The original renderer and its padded-root failure remain unchanged. No new
derivative/root probe, original scientific function, producer, eigensolver or
sampler is executed. This independent diagnostic uses Python standard-library
exact rational arithmetic locally; it is explicitly not a replay in the
original pinned Python/NumPy/SciPy runtime.

Embed the saved binary64 poles, weights, broadening, failed endpoints and
observed derivative padding as exact rational numbers. Check that every pole
is outside the failed bracket and every squared distance exceeds eta^2/3.
Consequently every Lorentzian curvature term is positive throughout it.
Let dmin/dmax bound each absolute pole distance. Using only pi>3 gives

M = (2 eta/3) sum w (3 dmax^2 - eta^2)/(dmin^2 + eta^2)^3,
B = (2 eta/3) sum w dmax/(dmin^2 + eta^2)^2.

M bounds the unrounded curvature from above throughout the retained bracket;
B bounds the sum of absolute unrounded derivative terms there. The exact
rational operations in this diagnostic require no numerical interval library.

Separately declare a **conditional** floating-evaluation model for the original
point enclosures: normal binary64 elementary arithmetic and small integer
powers, nearest binary64 pi, no underflow/overflow, and an accurate math.fsum
with at most two final rounding errors. The positive per-term path has fewer
than64 error factors (including subtraction and repeated denominator factors),
and fewer than1024 factors bound the complete11-term evaluation. A deliberately
loose absolute center-evaluation bound is rho=gamma1024*B, with unit roundoff
2^-53. This is a declared arithmetic assumption, not a hardware-verified
interval-arithmetic certificate. Checks retain wide normal-range margins.

If two endpoints have opposite *padded* signs with observed padding eps, their
computed derivative centers differ by more than2eps. The mean value theorem
and the conditional center bound imply that difference is at most
M times bracket width +2rho. Therefore, if M*tol+2rho <2eps, no subbracket of
width at most the original tolerance can have opposite padded signs within
the retained failed bracket under those assumptions. Merely trying more
quarter/midpoint probes would not satisfy the unchanged criterion there.

The conclusion is restricted to this saved bracket and arithmetic model. It
does not identify every root, alter the complete census, prove global solver
failure or establish that a different precision protocol will pass. Both the
original reference failure and zero pilot credit must remain visible.
