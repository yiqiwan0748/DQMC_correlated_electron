# Independent PySCF comparison: source and scientific evidence

This packet concerns exactly two catalogue questions: **P23**, two-orbital many-electron excitation spectra and spin labels, and **P24**, directed orbital-transition strengths. Seven finite fixtures, parameter values, individual numerical checks and serialization repairs are not additional benchmark ideas.

The native target calls the existing two-orbital kinetic/Hamiltonian APIs, with explicitly new dense-eigensystem, spin and directed-transition adapters. The independent reference calls actual PySCF 2.14.0 FCI Hamiltonian contractions, spin contractions and transition density matrices. This checks the retained native finite-cluster operators and complete low-energy manifolds; it does not test an iterative PySCF eigensolver or native DQMC dynamics, and it establishes no thermodynamic result. Only intraorbital Hubbard interaction is present.

## Retained failure and analysis repair

The original pilot job 6465868 succeeded in both scientific producers and then failed while serializing a NumPy Boolean in the comparator. Its original execution receipt remains **UNRESOLVED**, with process exits **[0, 0, 1]**, alongside the complete target/reference arrays and traceback. The first account-related sbatch rejection created no job; only a minimal derived history and the original receipt hashes are retained.

Job 6465910 reran analysis only, using the unchanged original comparator through a strict `np.bool_` to `bool` JSON adapter. The wrapper rejects other unsupported types. Original raw data, scientific formulas, thresholds and fixtures remain unchanged. The resource assessment retains the failed comparator cost, adds repair cost, and uses the original 2100-second/2800-MiB eligibility criteria. The independently reproduced forecast is 712.615451 seconds and 2641.671875 MiB. Neither pilot nor repair receives catalogue credit.

The [portable pilot replay](replay_pilot.py) verifies all three source archives and scientific data members, reruns the original comparator through the frozen strict serializer, and executes the unchanged scientific/resource fragment of the original independent auditor. It passes 122 selected scientific/resource checks. The original full independent audit passed 273 checks and remains identified by its hash; remote inode/transport and root-admission checks are deliberately not claimed by the public replay.

```bash
python -B replay_pilot.py --output-dir /a/new/empty/replay-directory
```

Use Python 3.11 or later and NumPy 2.4.6. Replay only reads retained eigensystems/operators; it does not invoke Julia, PySCF producers, scientific sampling, SSH or a scheduler.

## Frozen releases and publication selection

The [source index](sources/index.json) pins the original v1, analysis-only repair v2 and production v2 source archives and freezes. Production v2 changes only the comparator's NumPy Boolean conversion; its native producer, independent oracle, model, thresholds and resource policy are byte-identical to v1. The archives include the native source pins, official PySCF API references and licence. Installed binary hashes are metadata; no installed binaries are redistributed.

Original authored scientific `execution.json` and analysis process receipts are retained byte-for-byte. Root-admission bodies, whole prior queues, unrelated scheduler names, credentials, CV/contact information and transport/controller chatter are excluded. Per-dataset member manifests disclose the exact selected scientific files and original provenance hashes. New selection/replay metadata is explicitly labelled.

## Full production and independent public replay

Job 6465923 completed all seven registered fixtures and all 56 complete low-energy manifolds. All 3,503 frozen comparator gates pass. Its original independent audit passes 887 checks; a second replay directly from this compressed public payload passes 847 scientific/public-integrity checks with identical independently reconstructed case metrics. The difference in counts reflects the deliberately excluded private remote-inode/transport checks. Both audits and their exact source are retained.

```bash
python -B replay_production.py --output-dir /a/new/empty/full-replay-directory
```

This reruns only the frozen comparator and independent Fock-algebra audit on retained operators/eigensystems. It invokes no scientific producer or installed PySCF/Julia solver. `production/public_scientific_input_manifest.json` describes the role layout constructed during the recorded public replay; `public_frozen_comparator_replay.json` supplies its recorded `independent_replay/comparison.json` role. Fresh replay checks numerical fields with the original 1e-12 descriptive binary64 replay allowance, while retaining all original scientific thresholds.

Only the full-production outcome permits the conditional leaf builder to update P23/P24. A successful pilot, JSON repair or resource forecast alone never grants credit. Earlier v12 and its X01/P21 unresolved histories remain immutable.
