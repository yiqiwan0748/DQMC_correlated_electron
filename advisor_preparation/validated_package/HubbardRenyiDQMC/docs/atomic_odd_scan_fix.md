# Local proposal laziness repairs an atomic sampling failure

Date: October22,2026. This active package differs from the preserved source snapshot. The historical archive and accepted raw results were not changed.

At zero hopping, with an odd number of imaginary-time slices, the original deterministic time scan can have closed auxiliary-field orbits. For the five-slice pattern `(+,-,+,-,+)`, flipping each slice in order keeps `|q|=|sum_l s_l|=1` throughout the sweep. The per-site determinant weight is

`W(q) = 2 exp(beta*mu) [cosh(beta*mu)+cosh(alpha*q)]`,

where `alpha=acosh(exp(dtau*U/2))`. Since `W(1)=W(-1)`, every flip accepts. The reverse sweep retraces this orbit. More samples cannot correct its initialization-dependent class weight.

The failure was exposed by the external atomic benchmark `L=4,t=0,U=4,beta=0.5,dtau=0.1,mu_centered=0.3`. Seeds4001 and5003 initialized one site on the alternating orbit; the other six registered seeds initialized none. Exact enumeration of all32 fields at one site predicts total energy `-1.50960883266327` for one frozen site and `-1.88732299166887` when no site occupies that closed class. The two groups of saved chain means match these predictions. The correct unrestricted atomic energy is `-1.83079392864864`.

The complete forward-plus-reverse transition matrix has communicating-class sizes `[30,1,1]`. Its exact target is stationary, so checking stationarity alone would miss this failure. Uniform random initialization assigns mass `2/32=0.0625` to the closed orbit, whereas its correct equilibrium mass is approximately`0.037415239588219`. Within-chain blocking also misses the between-class problem.

## The correction

Before **every local HS proposal**, independently draw an idle decision with probability`1/2`. If idle, leave all state arrays, signs, determinants and tracker caches unchanged and proceed to the next scheduled site. Otherwise execute the existing proposal and Metropolis acceptance logic without modification. This changes each local kernel `P_l` to

`P_l_lazy = (I+P_l)/2`.

A convex combination of identity and a target-invariant local Metropolis kernel preserves the same target. Each local position now has a positive self-loop. In the finite positive-weight atomic fixture, a full scan can realize any individual bit flip by idling at every other position, making the corrected scan irreducible and aperiodic. The independently enumerated corrected32-state transition matrix has a single communicating class, positive diagonal probabilities and target-stationarity error below`1e-17`. Adding an idle decision only to an entire sweep would preserve the closed classes and would not fix this problem.

The change applies to exactly three active proposal loops:

- `_local_update_replica!`: ordinary scalar measurements and lambda0.
- `_local_update_replica_fused!`: direct nonzero-lambda replica weighting.
- `_local_update_replica_rank1_fused!`: optimized nonzero-lambda replica weighting used by the active node runner.

The public wrappers and all direct/optimized proposal formulas remain the same. The active node runner uses the corrected optimized path and delegates lambda0 to the corrected ordinary path. Instrumented historical timing scripts under `analysis/` and the immutable archive are provenance, not corrected production samplers.

## Statistics and compatibility

The returned `acceptance_rate` remains `accepted/N` over scheduled spatial slots. Idle slots contribute zero, so it now measures accepted moves per scheduled slot, not acceptance conditional on proposing a flip. It will usually be lower; no extra accepted-move claim is implied. Work per scheduled sweep is also lower because idle slots skip determinant/update arithmetic, while decorrelation per scheduled sweep changes. Runtime and effective-sample comparisons require fresh measurements.

The new idle draws change seeded trajectories. Direct and optimized active paths consume the same added idle draws and retain their same-seed equivalence test. Historical bitwise matching with the unmodified archived sampler is no longer expected. Historical performance numbers remain measurements of the historical implementation; they must not be relabeled as timings of this correction.

## Verification and reproduction

The standard package suite now includes `test/atomic_lazy_scan_tests.jl`. It checks the exact32-state transition graph, invariant target, predicted anomalous seed classes and correct atomic totals. It also exercises all three actual active APIs at `U=0`, where the original scan necessarily visits only two HS configurations, and verifies that the corrected APIs visit more than two while preserving exact free Green functions.

Run the targeted regression after preparing the package environment:

```powershell
julia --project=. --startup-file=no .\test\atomic_lazy_scan_tests.jl
```

The analytical proof uses only standard libraries and can run without package installation:

```powershell
julia --startup-file=no .\analysis\atomic_odd_sweep_proof.jl
```

`--kernel` additionally records exact-field Green-function errors, fresh stabilized reconstructions, reported stabilization corrections, propagator/field consistency, occupation bounds, alternating-orbit visits and trajectory fingerprints for stabilization periods1 and2. Its output path must be new. The upstream `JDQMCFramework1.2.6` constructor explicitly supports partial final stabilization intervals through `ceil(Ltau/n_stab)`; a five-slice grid with period2 is not an unsupported parameter combination.

The patched scientific benchmark arms must be rerun against the independent references. This code change alone does not certify all physical observables, convergence regimes or sign-problem cases.
