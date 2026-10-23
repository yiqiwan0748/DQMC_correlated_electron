# Prospective Monte Carlo execution plan

Registered before the first interacting DQMC arm starts on ISAAC, 22 October 2026.

The 24 stochastic comparisons are the 16 QT and eight ATOM cases in the canonical
manifest. Numerical tolerances and equivalence margins remain fixed. The eight
seeds, 512 burn-in sweeps, 64 bins, four lambda increments and estimator are fixed
in the reviewed consumer.

First time one QT-01 scalar arm at 512 measurements per chain. This is a runtime
pilot; its result will remain separate from the final campaign. Then run every
stochastic comparison at 4,096 measurements per chain. If a comparison remains
unresolved through inadequate precision or blocking, the only prespecified
extensions are 16,384 and 65,536 measurements per chain. Retain every earlier
receipt. A numerical discrepancy is investigated, not hidden by repeating seeds.

Each interval uses two-sided Student t7 alpha=0.01/100. There are at most
24*3+1=73 stochastic interval inspections in this plan, so a union bound remains
below 1%, under the independent-chain normal approximation stated in the guide.
The same seeds or dependence between inspections do not invalidate that union
bound. Diagnostics do not prove equilibrium or normality. Failed blocking or
effective-bin checks remain unresolved.

The 65,536-measurement stage is a precision limit, not a guarantee of scientific
completion. Other models, larger systems, long-time instabilities and continuum
or thermodynamic limits require their own external benchmarks.
