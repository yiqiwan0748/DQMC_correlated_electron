# Atomic-limit sampling failure and repair review

Status: **the original deterministic local-proposal schedule fails an exact atomic-limit benchmark; all 24 corrected atomic/thermal external comparisons now pass their registered criteria**. The first corrected run and its unresolved precision cases remain visible below. These finite-fixture results are recorded separately from the analytical repair proof and do not certify every physical regime.

The new external benchmark exposed a sampling defect that ordinary within-chain blocking did not detect. Preserve the original source, original results and this diagnosis together. A historical numerical agreement is evidence for its stated calculation, not a certificate that every model limit is sampled correctly.

## Observed ISAAC result

Original receipt: [external_atomic_mc_4096.toml](../isaac/results/external_atomic_mc_4096.toml), job **6456805**, 22 October 2026, eight independent seeds, 512 burn-in cycles and 4,096 measured cycles per seed, retained as 64 bins. Its SHA-256 is `a69af55fb680a1b44eec95cdbae8bdb77439a456e9352031c660ca7ab323e094`.

The affected fixture is four independent Hubbard atoms: `L=4,t=0,U=4,β=0.5,μ_centered=0.3,Δτ=0.1`, with **five** time slices. The centered physical Hamiltonian is `H=U Σ(n↑−1/2)(n↓−1/2)`; `−μ_centered N` belongs to the sampling weight and is excluded from the reported physical energy. Thermal ED and the four-state atomic formula agree. The original DQMC energy, density and double occupancy remained unresolved; the entropy arm's old precision gate passed. That entropy label does not certify this affected kernel.

| Quantity | Exact four-atom target | Original DQMC mean |
|---|---:|---:|
| Centered physical energy | −1.830793928649 | −1.801098546593 |
| Total occupancy | 4.161481159435 | 4.163691758412 |
| Total double occupancy | 0.623042097555 | 0.631571242558 |

Seeds `4001` and `5003` produce energy means `−1.504941` and `−1.512177`; the other six lie between `−1.906280` and `−1.889467`. Each chain stays near its own plateau throughout the saved bins. The between-chain energy standard error is `0.0638783`, versus a naive pooled within-bin estimate `0.00385368`. Classic unsplit Rhat on the saved energy bin means is **2.297285**; a control fixture with ten time slices gives **1.000743**. These are diagnostic statistics, not confidence probabilities. The old estimated 469 effective bins does not detect the inaccessible classes.

[Read-only bin analysis](atomic_bin_diagnostic.json) records per-seed means, early/late halves, blocking, autocorrelation of saved bins and exact measurement identities. Sweep-level autocorrelation cannot be recovered from bins averaging 64 cycles. Atomic energy/density identities hold to `6.5×10⁻¹⁵`; off-diagonal one-body terms vanish to numerical precision. Consistent measurements do not establish correct sampling.

## Exact cause

For one isolated atom, let `s_l=±1`, `M=β/Δτ=5`, `λ=acosh(exp(Δτ U/2))` and `x=exp(β μ_centered)`. After the discrete spin Hubbard–Stratonovich transformation, the field weight, up to a state-independent factor, is

`w(s)=[1+x exp(λ Σ_l s_l)] [1+x exp(−λ Σ_l s_l)]`.

The two alternating temporal patterns `(+−+−+)` and `(−+−+−)` have `|Σs|=1`. A deterministic forward scan proposes every field flip once; every intermediate sum still has magnitude one. Each Metropolis ratio is one, so every flip is accepted. One forward scan exchanges the two patterns; a following reverse scan returns to the original pattern. The package's measured forward-plus-reverse cycle therefore has **two separate fixed states**, not one communicating two-state cycle.

Independent enumeration of all 32 single-atom field states gives:

| Transition convention | Original communicating-component sizes | With per-local proposal probability 1/2 |
|---|---|---|
| One forward scan | 30, 2 | One component of size 32 |
| Forward then reverse cycle | 30, 1, 1 | One component of size 32 |

The closed patterns have combined target probability `0.0374152395882`, whereas a uniform initial field assigns them probability `2/32=0.0625`. Because the scan cannot enter or leave those classes, more cycles cannot repair this class-population error. The four-atom means conditional on no trapped site are `E=−1.887322991669,N=4.157272993710,D=0.606805748938`; with one trapped site they are `E=−1.509608832663,N=4.185390986342,D=0.715293285005`. These independently derived plateaus explain the observed seed split. [Seed initialization analysis](../isaac/diagnostics/atomic_odd_sweep_diagnostic.toml) separately identifies an alternating site for seeds 4001 and 5003; the original saved receipts themselves do not store full trajectories.

The corresponding infinite-chain-length, uniform-initialization class mixture would have `E=−1.792894451917`, rather than the exact `−1.830793928649`. This is an analytical expectation over initializations, not a measured run result.

[Exact enumeration code](exact_atomic_scan_graph.py) and [machine-readable proof](exact_atomic_scan_graph.json) import no user sampler. Reproduce with `python3 external_oracles/exact_atomic_scan_graph.py` from the `advisor_preparation` directory. For the numerical fixture, both original and repaired transition matrices preserve the target distribution to `1.4×10⁻¹⁷`; stationarity alone is insufficient because the original matrix is reducible. The independent source is the [discrete Hubbard–Stratonovich transformation](https://doi.org/10.1103/PhysRevB.28.4059); the finite transition-graph calculation and conditional means are derived here.

This diagnosis establishes a defect at `t=0` with the stated odd-time-slice schedule. The alternating construction also applies to isolated atoms with odd `M≥3`, positive repulsion and balanced spin chemical potentials. It does not establish biased historical results at nonzero hopping, nor a general odd/even classification for all interacting lattices. A stabilization interval not dividing `M` is not the cause demonstrated here.

## Mathematical review of the correction

The active correction inserts an independent idle probability `1/2` **before each local proposal**, in the ordinary, direct fused and optimized rank-one fused paths. Imaginary-time propagation and stabilization remain outside this local guard. The review artifacts are [the source patch](../isaac/diagnostics/atomic_local_lazy_fix.patch) and [source hashes](../isaac/diagnostics/atomic_local_lazy_fix_receipt.json). Acceptance statistics count idle slots in their denominator; their values therefore have a changed interpretation.

For strictly positive target weight `π(s)`, a proposed flip from `s` to `s'` now has transition probability

`P(s,s') = p min[1,π(s')/π(s)],  p=1/2`.

Detailed balance holds for each local kernel because `π(s)P(s,s')=p min[π(s),π(s')]`. The remaining probability is a hold. Composing these local kernels in a fixed scan preserves `π`; that composition need not itself satisfy detailed balance. With all field configurations of positive weight, one full scan can flip exactly any one selected field and hold every other proposal with positive probability. The full-scan chain is therefore irreducible on the finite hypercube and has a positive self-loop, so it is aperiodic. The same conclusion holds for a forward/reverse cycle. Exact enumeration confirms one 32-state component after the correction; target-stationarity residuals are below `7×10⁻¹⁸`.

The proof assumes connected strictly positive field support. It does not certify ergodicity in every determinant-zero or sign-changing many-body problem. A single idle step at the **whole-sweep** level would retain the same disconnected classes. Permuting only spatial sites cannot release a trapped isolated atom. The new local RNG draws intentionally change trajectories; historical same-seed comparisons and speedups apply to their historical source until freshly remeasured.

## Required evidence after the correction

1. Execute the actual generic, direct fused and rank-one fused paths with forced alternating initial fields, including an atomic odd-slice fixture, and check fresh Green matrices. Exact transition enumeration is a separate proof, not a substitute for invoking those paths.
2. Rerun the externally anchored atomic observables and finite-hopping thermal comparisons from fresh initial states. Retain all seeds, bins, source hashes, controls, timings and failed/unresolved outcomes. Do not combine old-kernel PASS rows with a corrected-kernel summary.
3. Apply the preregistered precision criteria and an explicit between-chain mixing gate. [Retrospective diagnostics](retrospective_split_rhat.json) recompute classic split-Rhat from paired jackknife pseudovalues without changing any original receipt. A diagnostic pass is not proof of equilibration. Modern [rank-normalized and folded Rhat](https://arxiv.org/abs/1903.08008) can diagnose failures that classic Rhat misses; it is a useful further extension.
4. Treat Student-t intervals across eight independent chains and associated multiple-inspection allocation as nominal, conditional on sufficiently equilibrated chain estimates and a suitable sampling approximation. Arithmetic alpha allocation does not prove finite-run Markov-chain coverage.

No thresholds should be relaxed to turn the original unresolved cases into passes. Keep the original atomic failure visible after a successful corrected run, because it explains the source change and demonstrates why external model limits and independent chains are necessary.

## Recorded supporting-test execution

ISAAC job `6456921` passed 382 supporting checks: 254 package workflow checks, 14 exact atomic-transition checks and 114 actual generic/direct-fused/rank-one-fused free-field checks. It then stopped on a Julia world-age error in the additional diagnostic harness, before running the external atomic comparisons. The [failed log](../isaac/logs/atomic_corrected_6456921.log) is retained. These checks support the implementation change; they do not substitute for the 24 corrected external stochastic comparisons. The [ISAAC reproduction guide](../isaac/ISAAC_REPRODUCE.md) preserves the original failed wrapper and uses a versioned diagnostic/resume wrapper without replacing the physics source archive.

## First corrected atomic execution

Job `6457036` used the versioned v2 diagnostic, then ran the eight external atomic comparisons with the corrected frozen kernel. The [field/Green diagnostic](../isaac/results/atomic_corrected_kernel_diagnostic_v2.toml) records 16 seed/stabilization arms. Every site revisited the alternating patterns, with at least 30 visits over the retained diagnostic measurements. Every stabilization-period-two trajectory matched its period-one counterpart. Maximum equal-time Green error against the exact field formula was `5.56×10⁻¹⁶`; comparison with fresh stabilized reconstructions gave zero at the checked cuts. These are observed diagnostic results, not an extrapolation to every lattice.

At 4,096 measurements per chain, the [first corrected external atomic receipt](../isaac/results/external_atomic_corrected_4096.toml) has **7 PASS, 0 FAIL, 1 UNRESOLVED**. All eight split-Rhat values meet 1.05. `ATOM-05` is precision-limited: absolute energy discrepancy `0.01021856` plus simultaneous half-width `0.02492871` exceeds its unchanged `0.03` margin. Its new split-Rhat is 1.0, compared with 2.22692 for the original trapped result. Longer sampling remains subject to the registered extension plan; this initial result is not overwritten.

The [independent receipt audit](corrected_atomic_independent_audit.json) checks the kernel/consumer hashes against the frozen corrected archive and recomputes all eight mean/error/Rhat acceptance decisions from the saved data. The audit itself launches no scientific jobs. Passing seven finite atomic targets does not certify unsupported observables, arbitrary determinant-zero/sign-changing supports, low-temperature stability or continuum/thermodynamic limits.

## Final corrected finite-fixture outcome

The registered precision extensions resolve `ATOM-05` and `QT-09,QT-15,QT-16` at 16,384 measurements per chain. The [atomic extension](../isaac/results/external_atomic_corrected_16384.toml) and [thermal extension](../isaac/results/external_thermal_corrected_16384.toml) are separate receipts; the initial 4,096-measurement results remain unchanged. Combining the latest corrected result for each fixture gives **8/8 atomic and 16/16 thermal DQMC comparisons passing**, with all original margins and the split-Rhat/blocking/effective-bin gates retained.

The [final independent audit](final_campaign_independent_audit.json) reconstructs paired jackknife chain means, standard errors, Student-t7 half-widths, split-Rhat, blocking changes and influence-bin effective counts from all 28 corrected inspections. It checks kernel/consumer hashes against the frozen corrected archive. Combining these 24 latest stochastic arms with 164 unchanged deterministic solver arms independently reproduces **100 passing finite-system comparison cases and 188 passing method arms**, matching the [campaign summary](../isaac/results/summary/summary.json). No original-kernel stochastic PASS enters this corrected result.

This establishes the tested atomic and small matched thermal observables, alongside the separately recorded finite ED/DMRG/free-fermion comparisons. It does not establish general ergodicity, thermodynamic-limit accuracy, universal error-bar coverage or reliability for unsupported models. The original reducibility proof, failed chains, first corrected unresolved cases and harness failure remain part of the record.
