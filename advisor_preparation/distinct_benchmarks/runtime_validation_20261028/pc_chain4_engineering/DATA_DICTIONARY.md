# Data dictionary

All tables select values already stored by the completed run. Decimal text round-trips to the original binary64 numbers, including signed zero. Site indices and prefix cuts are one-based. Complex values retain separate real and imaginary components. A prefix cut b separates sites 1 through b from the remaining sites. S2 uses the natural logarithm, S2 = -log Tr(rho_A squared), for the normalized pure state.

| File | Rows | Meaning |
| --- | ---: | --- |
| `data/stage_scalars.csv` | 3 | Snapshot raw norm squared; native energy numerator, explicit norm denominator and their stored ratio. These are same-state measurements. |
| `data/entropy.csv` | 9 | Native user-helper S2 and its norm check at three cuts per snapshot; separately captured SVD purity numerator, denominator, ratio and adapter S2. Keep both floating observations; do not silently equate them. |
| `data/spin_correlations.csv` | 48 | All 4 by 4 spin-dot entries at each snapshot. Contains the separately measured raw numerator and normalized ratio plus native correlation-API SzSz, S+S-, S-S+ and spin-dot entries. Matrix indexing is column-major in the raw JSON. |
| `data/local_observables.csv` | 60 | Nup, Ndn, Ntot, Nupdn and Sz for all four sites and three snapshots, with stored numerators, denominators and ratios. |
| `data/optimizer_ritz_energies.csv` | 12 | Final-bond optimizer Ritz diagnostic for each completed sweep. This is not the post-truncation saved-state energy. |
| `data/bond_diagnostics.csv` | 72 | Both directions of each sweep: bond, retained dimension, signed native truncation value, actual controls and local optimizer energy. No absolute-value conversion or replacement was applied. |
| `data/stop_callbacks.csv` | 12 | The original stopping callback decisions. |
| `data/scientific_raw.json` | Structured | A scientific projection of the saved native-arm JSON. Contains all retained arrays, tensor/state witnesses and observations; removes private execution metadata and unused opaque index text. Numeric values are unchanged. |
| `data/recorded_exact_replay.jsonl` | 2,838 | Byte-identical scientific replay log: 2,799 checks, 17 state records, 21 entropy-cut records and one completion record. This includes exact rational diagnostics and entropy enclosures; no reference replay was repeated to make this file. |

The spin convention is S_i dot S_j = Sz_i Sz_j + (S+_i S-_j + S-_i S+_j)/2. The local occupation order is empty, up, down, up-down. The independent fermion reference uses all up modes first, followed by all down modes, with occupations and signs explicitly retained. Tensor amplitudes use site 1 as the fastest local index.

The observer's zero-noise window at sweep 3 is explicitly `NOT_APPLICABLE_NO_ZERO_NOISE_SWEEPS` with a null value, not a zero estimate. Later zero-noise windows and every raw observer maximum are retained in the scientific JSON. A callback maximum and a saved-state error are different quantities; no ground-state error is inferred from either.

`code/compute_s2_itensor_dmrg.jl` is the unchanged general-purpose original analysis program. Its defaults are not a reproduction command for this engineering arm. `code/EngineeringCompare.jl` is the unchanged exact-reference module used by the completed run. It expects the documented raw structure through a separately supplied loader/callback; this packet does not supply a tested standalone runner. Both files are source for inspection, and neither was executed during packet preparation.
