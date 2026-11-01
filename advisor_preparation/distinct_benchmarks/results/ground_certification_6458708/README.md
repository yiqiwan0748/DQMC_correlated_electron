# Initial ground-certification attempt6458708

Independent result: **X12 PASS; X09 UNRESOLVED**. This attempt adds one scientific completion. A successful audit of this mixed outcome is not a claim that both benchmarks passed. All original raw receipts and the frozen version1 protocol remain unchanged.

| Native DMRG start | Energy above independent ground | Positive residual-norm variance | Ground-projector fidelity | Certificate |
|---|---:|---:|---:|---|
| neel |2.81755e−8 |1.85606e−8 |0.999999909502 | Unresolved |
| charge_density_wave |4.05063e−8 |2.84782e−8 |0.999999870372 | Unresolved |
| spin_domains |2.31548e−9 |2.31617e−9 |0.999999992825 | Pass |
| shuffled |5.06945e−8 |3.55633e−8 |0.999999838151 | Unresolved |

The unchanged limits are energy error1e−8, variance at most1e−8 with a−1e−9 roundoff allowance, and fidelity at least0.99999999. Three starts miss these gates. Their variance was independently recomputed as ||(H−E)psi||²/||psi||², establishing a positive residual rather than moment-subtraction roundoff. All states have the correct norm and particle sector. Reported truncation errors become zero after sweep3, but that alone does not certify the ground state.

The source ED Hamiltonian equals the independent QuSpin matrix exactly; QN and dense user MPO matrices agree within6.6614e−16. Direct state-amplitude Rayleigh moments and fidelities reproduce the frozen comparator. The manufactured first-excited and highest-eigenvector MPS controls have residual variances about1.99e−29 and2.47e−28 and ground fidelities around1e−18; both correctly fail ground certification. They are deliberately constructed controls, not observations of natural optimizer trapping.

This evidence identifies insufficient precision in this specified eight-sweep run, rather than a demonstrated Hamiltonian or package bug. The user production command's default is14 sweeps; the initial benchmark adapter used8 and the native library's default local Krylov dimension3/rebuild limit1. Any higher-effort retest must be a new declared release with these receipts preserved and every scientific threshold unchanged.

X12 independently passes all ten pure/mixed/rotated/canonical-Gibbs prescriptions and both subsystem cuts in each backend. The largest user-versus-analytic RDM-spectrum L1 error is9.0975e−16 and purity error1.3323e−15. Pure/mixed one-site S2 differs by0.28768207 and connected Czz by0.0625, confirming that the fixture distinguishes the ensemble prescriptions. Complex transverse-spin means, complete canonical-N3 states, full ground-mixture cooling and basis-rotation invariance are verified.

The fresh frozen comparator replay is byte-identical. `independent_audit_attempt01.json` binds23 inputs, the eight-file source archive, all five native source/environment pins and six remote/local text transfers. Its exact scalar/vector auditor is copied to `independent_auditor_attempt01.py`. Scheduler result: COMPLETED/0:0,5m53s. No new model, eigensolver, optimizer or Monte Carlo was run during this audit.
