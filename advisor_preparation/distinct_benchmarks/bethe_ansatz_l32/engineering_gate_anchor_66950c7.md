# Independent Bethe pilot audit and next resource bridge

ISAAC job 6458253 completed successfully in 6m03s. All 18 downloaded scientific receipts/log files match their remote SHA256 values; all 12 separate solve receipts equal their aggregate target rows. The immutable pilot source closure is unchanged. The frozen summary replay matches exactly. The independent audit passed 149 provenance, normalization, sector, variance, truncation, symmetry, quadrature-accounting, and convergence checks.

**A18, A19, and A20 remain NOT_RUN in the scientific catalogue.** The raw pilot summary correctly labels them UNRESOLVED/PILOT_ONLY and adds zero scientific completions. No source file, existing result, canonical ledger, or job was modified to prepare this audit. Local activity was read-only inspection and arithmetic/postprocessing, not a new physics run.

| Pilot quantity | Measured value |
|---|---:|
| Oracle process wall time |21.82s|
| Julia target process wall time |338.80s|
| Target internal elapsed time |239.826s|
| Sum of 12 per-solve times |236.685s|
| First solve |131.202s|
| Remaining 11 solves combined |105.483s|
| L12 refined solve range |18.967–19.088s|
| L12 coarse solve range |8.814–10.029s|
| Target `/usr/bin/time` peak RSS |1,307,856 KiB, about 1.25 GiB|
| Slurm accounting MaxRSS |1234M|

The two memory figures are preserved as recorded; Slurm polling and the process resource monitor need not capture identical peaks. First-use latency and process/IO overhead are substantial. JIT compilation and filesystem wait were not separately instrumented, so the entire first-solve delay cannot be assigned to either. Do not scale total RSS by a power of bond dimension or treat a warm small-chain solve as a cold-job estimate.

The maximum coarse/refined differences are `1.025e-9` in energy per site, `8.711e-8` in the charge-gap combination, and `1.152e-9` in mean double occupancy. The maximum particle-hole energy error is `2.665e-15`. L12 refined variance per site is at most about `6.0e-11`, and final truncation error at most `4.63e-12`. Both L12 arms reach their bond caps 128 and 256. This demonstrates the pilot's convergence at these finite sizes, not an error bound for larger systems.

Finite-size effects remain substantial. For the refined L12 chain, the energy/site differs from the bulk integral by `+0.029876`, the charge gap by `+0.425458`, and mean double occupancy by `-0.005458`. These are expected finite-chain offsets, not failed solver comparisons. Two pilot sizes do not support the predeclared controlled thermodynamic extrapolation. The original full proposal's training sizes and held-out sizes must remain unchanged by this planning audit.

## Recommended next bounded bridge: L24 first

Prepare a **new pilot-only release**, leaving the original pilot frozen. Reuse the same hashed user MPO and the same additive target logic; do not edit the kernel. Keep U4,t1,OBC and the three half/minus/plus sectors. At L24 use two independently seeded arms with the prospective full-execution settings: coarse chi 256, 18 sweeps, cutoff 1e-12; refined chi 512, 24 sweeps, cutoff 1e-13. This gives six solves and directly measures the bond cap that the proposed full programme would need. A same-cap128/256 bridge would be cheaper, but would leave chi512 memory and runtime unmeasured.

Use one CPU and a hard 3500M memory limit. A conservative preliminary engineering range is **15–40 minutes after start**, with a **45-minute wall cap** and a 40-minute target-command cap. This is not a measured L24 runtime. A crude `L*chi^3*sweeps` estimate applied to the measured 19-second L12 refined solve suggests about 405 seconds for one L24 refined solve at doubled cap and 24 sweeps; QN block sizes, actual retained bonds, contractions, variance evaluation, cold start and filesystem latency can change that markedly. The range and hard cap are therefore deliberately broader than that crude extrapolation. Scheduler wait is separate.

Retain the original pilot's normalization, sector, symmetry, variance, truncation and observable-stability tolerances. Do not loosen a gate to make a bridge finish. Record every solve's norm, Nup/Ndown, contracted energy, variance, complete truncation history, site double occupancy, final/max bond dimension and elapsed time; retain all completed per-arm receipts if the job hits a cap. Continue reporting zero scientific completions. A memory limit, timeout, or convergence failure is a measured feasibility limitation and requires a separately reviewed new resource/precision protocol.

Before root submits, freeze the exact new manifest/source archive, check the live queue, and state an EDT check point. At most two scientific jobs may run in total, including unrelated jobs; keep this bridge in the serial secondary lane. No bridge submission or executable-source alteration is part of this audit.

## L32 follows only after the L24 receipt review

Do not queue L32 automatically now. It should use the same three sectors and intended chi256/512 arms only if all six L24 cases pass the unchanged quality gates and the resource evidence leaves room. A conservative prerequisite is L24 job wall time at most 30 minutes and peak RSS at most 2300 MiB; these are engineering continuation criteria, not scientific tolerances. If either criterion fails, revisit allocation/implementation before a larger pilot.

For an eligible L24 result, a proposed L32 bound is one CPU, 3500M and a 60-minute wall cap. Its runtime range must be recalculated from actual L24 timings rather than asserted now. The size ratio 4/3 at fixed chi is only a planning factor; entanglement, effective QN blocks and convergence may dominate it. The hard memory cap remains authoritative even if a linear size estimate looks favorable.

The L24/L32 bridge checks resource and truncation feasibility. It neither replaces the held-out L56/L64 checks nor creates a bulk result. Before any full scientific execution, root must separately implement and freeze the full extrapolation/uncertainty procedure from the predeclared proposal, including all training sizes, model/drop-size checks, held-out prediction gates and propagated sector-energy uncertainty. Bridge data must not be used to change those scientific choices. Any reuse of prior data must be declared in that new release before full analysis; the simplest policy is to keep resource-bridge receipts as separate diagnostics and obtain the required full-release arms explicitly.

## ALF status and retrieved evidence

Job 6458258 completed in 1m11s, exit 0. Actual status files say `ALF_BUILT_NOT_BENCHMARKED` and `SETUP_BUILT_NOT_BENCHMARKED`; the compiler log ends in successful links. Appended `-fno-fast-math` overrides the earlier default `-ffast-math`. This is setup evidence only: P21 remains NOT_RUN. Text/JSON/log receipts are saved separately under `results/alf_build_6458258`; no binary, system library, compatibility symlink, or complete ALF source checkout was downloaded.
