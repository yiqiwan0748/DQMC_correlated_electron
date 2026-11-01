The original P16 comparison failed for the disjoint spatial set `{0,2}` on the four-site, periodic, half-filled Hubbard chain at `U/t=4`. The actual user ED state has energy −2.1027484834620633 and residual 2.04×10⁻¹⁴. The independent QuSpin state has the same energy within 1.8×10⁻¹⁵ and residual 1.65×10⁻¹⁴. Single-site purities agree within 3.4×10⁻¹⁶. The union purity differs: user 0.18328212068756122, QuSpin 0.19215599964073315. The resulting Rényi-2 mutual-information discrepancy is 0.047280931427645. The immutable failure is [packages.json](../../results/tranche01_attempt01/packages.json).

The existing `renyi2_entropy_sites` applies the sign `(-1)^(Nup_B*Ndown_A)` when regrouping the fixed-sector Fock basis into A and B. For arbitrary spatial sets, the full fermionic regrouping also moves occupied B modes past occupied A modes *within each spin species*. Writing `Inv_sigma` for the number of occupied `(i in B,j in A)` pairs with `i<j`, the required sign is

```
(-1)^(Inv_up + Inv_down + Nup_B*Ndown_A).
```

This follows directly by counting inversions between the original mode order `up(0..L-1),down(0..L-1)` and `up(A),down(A),up(B),down(B)`, with each spatial list sorted. The omitted within-species parity generally depends jointly on A and B, and thus changes the singular values of the coefficient matrix for disjoint cuts.

Contiguous prefixes have both omitted inversion counts zero. For a contiguous suffix in a sector with fixed Nup and Ndown, each extra count is `N_sigma_A*(N_sigma-N_sigma_A)`, a phase depending only on subsystem particle number, which is a local change of basis and preserves the entropy. This diagnosis therefore does not invalidate the prior prefix/suffix entropy results. It does identify an error in the claimed support for general spatial sets.

QuSpin's [spinful partial-trace implementation](https://quspin.github.io/QuSpin/modules/quspin/basis/basis_1d/fermion.html) converts the specified spin-up/spin-down spatial lists into the full fermion-mode list and calls its fermionic partial trace; its basis declares anticommutation for all modes. The diagnostic compares every four-site spatial subset to that independent implementation, rather than assuming a tensor-product trace with no fermionic shuffle.

`regrouping_v2.jl` is an isolated candidate function around the actual user's ED state. It leaves the frozen solver and failed receipts unchanged. `verify_quspin_regrouping.py` separately checks all 256 Fock configurations and 16 spatial sets against explicit occupied-mode inversion counting, then compares all 16 ground-state purities against QuSpin. `recheck_user_v2.jl` recomputes the actual user ED state, records old and candidate results, and compares P16 to the frozen external oracle at the original 10⁻⁸ tolerance. These are diagnostics and a correction attempt for P16, not extra benchmark ideas. Actual execution receipts are required before any corrected PASS claim.
