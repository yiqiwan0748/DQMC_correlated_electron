# Pinned PySCF reference environment

This setup installs the independently maintained PySCF 2.14.0 full-configuration-interaction reference and pinned dependencies into a new ISAAC environment. No Hubbard solver, model fixture or acceptance comparison runs in this setup, and it adds zero catalogue completions. The existing QuSpin and TeNPy environments are preserved.

Only version metadata, wheel hashes, installer source and execution receipts belong in the public record; third-party binary wheels and the runtime stay on ISAAC. New P23/P24 source/protocol releases and resource pilots must be reviewed and frozen before scientific execution.

Reference documentation: [PySCF FCI](https://pyscf.org/pyscf_api_docs/pyscf.fci.html), [PySCF source and release](https://github.com/pyscf/pyscf/tree/v2.14.0).
