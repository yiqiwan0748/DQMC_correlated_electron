# Corrected sampler campaign, registered before execution

22 October 2026, after the initial 4,096-measurement campaign exposed reducible
atomic-limit sampling. Original receipts and source remain immutable. The exact
32-state transition-graph investigation is a supporting correctness test, not an
additional external benchmark case.

The original systematic sweep has closed communicating classes at the second
atomic fixture. Longer runs of that kernel cannot establish reliability there.
The corrected working copy introduces an independent 50% idle probability at
**each local auxiliary-field proposal**, in ordinary and coupled-replica paths.
Each local kernel retains its target distribution; their composition retains the
target too. The full sweep need not itself satisfy detailed balance. The exact
atomic transition matrix must be irreducible, aperiodic, and stationary for the
analytical HS weights. This is a proof for that finite fixture, not for every
fermion model or every zero-weight surface.

Rerun **all 24 stochastic cases**, including previously passing ones, with the
corrected source, eight original seeds, 512 burn-in sweeps, 64 bins, and the same
four lambda increments. First use 4,096 measurements per seed. Fixed optional
precision extensions are 16,384 and then 65,536; select every unresolved case at
each stage. Never repeat a stage with new seeds or relax a physics margin.
Retain every earlier result. The corrected summary combines unchanged
deterministic solver arms with corrected Monte Carlo arms, never with old-kernel
Monte Carlo passes.

Add a split-chain Rhat gate of at most 1.05 on paired jackknife pseudovalues,
alongside the existing blocking and effective-bin gates. Scalar pseudovalues use
paired sign reweighting; entropy pseudovalues use paired leave-one-block-out
log-ratio estimates across nodes. This catches between-chain discrepancies that
the original within-chain diagnostics missed. Rhat is a diagnostic, not a proof
of equilibrium, independence, or approximate normality. Gate failure remains
UNRESOLVED even if the mean agrees with an oracle.

Corrected confidence intervals use Student t7 with alpha=0.01/400 per inspection.
The original campaign made 25 inspections (24 cases and one timing pilot), at
alpha=0.01/100. At most 72 corrected inspections are permitted. Their union-bound
nominal budget is 25*0.0001+72*0.000025=0.0043, below 1%, under the independent-chain
normal approximation already stated. Finite MCMC means need not follow an exact
Student distribution; this is not a rigorous global coverage guarantee. This
conservative accounting does not
repair bias from a reducible or unequilibrated original chain; those original
results are retained as diagnostic history.

Measured original runtimes were 4m41s for the 16 thermal comparisons and 1m42s
for the eight atomic comparisons, including compilation. Expect roughly 2–6
minutes per corrected 4,096-measurement job; the fixed extensions scale up to
about four and sixteen times their sampling portions. Queue delays are separate.
The 65,536 stage is a precision limit, not a promise of passing results.

The diagnostic rationale follows the multiple-chain variance comparison of
[Gelman and Rubin (1992)](https://dash.harvard.edu/entities/publication/73120378-8481-6bd4-e053-0100007fdf3b).
The particular paired pseudovalue implementation and threshold here are explicit
campaign choices, not claims of a general convergence certificate.
