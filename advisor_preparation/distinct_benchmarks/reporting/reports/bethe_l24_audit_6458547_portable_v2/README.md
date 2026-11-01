# Portable L24 feasibility audit, revision2

Six native user-MPO/ITensor DMRG solves pass the unchanged convergence gates. All12 original remote scientific/log hashes agree. A fresh replay of the exact summarizer extracted from the frozen source archive agrees byte-for-byte with both the original raw summary and preserved replay. This is arithmetic/postprocessing only; no numerical physics was rerun.

The original auditor and first report are preserved. Required original hash/scheduler and replay inputs are now named sidecars in the raw result directory; no temporary input path is required. The two historical input sidecars retain exactly their original bytes.

The30-minute/2300-MiB L32 engineering criteria are bound to the exact previously published README at commit66950c791f7f939ca176e5427be62cd8dfca8348 (2026-10-22T19:12:58Z), preserved as a raw sidecar with its Gitblob and SHA256. Job elapsed time16m31s=991s and target peakRSS1,484,548KiB=1449.75390625MiB satisfy those criteria. They are resource continuation gates, separate from scientific convergence tolerances.

All six sector/arm receipts and complete site/truncation arrays remain in the raw directory. Coarse/refined charge gaps differ by about1.87e-7; the refined finite-size gap1.45911 still differs substantially from the infinite-chain oracle. This permits a separately frozen larger feasibility pilot; it establishes no thermodynamic result. A18/A19/A20 remain NOT_RUN and scientific completions remain0.
