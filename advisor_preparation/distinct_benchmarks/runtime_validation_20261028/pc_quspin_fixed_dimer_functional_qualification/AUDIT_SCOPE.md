# Audit scope

The outcome is a functional qualification of one recorded QuSpin reference process for the fixed two-site model in the protocol. The saved audit checked the reviewed scientific source, complete output and field links, recorded scalar gates, and the required process qualification evidence. Publication preparation copied and checked saved data; it performed no new physical calculation.

The scientific routine registered twelve explicit constructor/method entries: one basis constructor, four number-operator constructors and four dense exports, then the Hamiltonian constructor, dense export, and complete eigensolver. The successful check retained all supported finite raw number matrices, Hamiltonian entries, four energies, and sixteen eigenvector entries.

The analytic reference was fixed before the data. An integer square-root enclosure and rational endpoints define the four reference levels. Each reported binary64 energy is represented exactly as a dyadic rational and checked against both endpoints of its reference interval. The limit is an absolute energy error of \(10^{-11}\). The stored-matrix Hermitian-defect limit is \(10^{-12}\). The publication checker validates the recorded scalar bounds; it does not generate a new reference.

The saved residual and Gram-defect matrices are point diagnostics of the stored binary64 data. They do not certify an ideal normalized eigenvector, ground-state selection, a global gap, entropy, correlations, or other observables. No native-solver comparison was performed.

These are finite floating-point outputs. The energy comparison is conditional on the registered model, operator convention, and analytic reference. It is not a general floating-point correctness proof for an operator, solver, state, or observable.

Two defects were found before numerical execution:

- The first source revision did not enforce every API witness predicate at both numerical boundaries. The corrected validation applies the same per-witness conditions at each required stage.
- The second source revision rejected legitimate empty Python module files. The correction accepts exact registered empty general module files while retaining positive-size requirements for API definitions and compiled images.

These were pre-use validation defects. The numerical check used the corrected conditions. This history supplies no additional numerical runs or scientific credit.

Malformed, unsupported, or nonfinite returned objects cause failure. The serialization contract covers supported valid arrays; an unsupported eigenvector object can prevent saving its paired returned energies. This successful record contains the complete supported finite arrays. Arbitrary foreign-object serialization is not claimed.

The public files contain scientific sources and results. Machine-specific provenance and execution envelopes are not included. The qualification does not establish a unique internal BLAS/LAPACK/JIT implementation, general package reliability, or production readiness.

The user's ED, DQMC, and DMRG implementations remain outside this test. Full X08 remains **NOT_RUN**. New ideas remain **0**, and the catalogue count is unchanged.
