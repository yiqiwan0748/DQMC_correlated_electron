# Initial external comparison batch

**100 cases passed; 188 solver arms passed; zero failed or unresolved arms.**

[Case table](case_status.csv) · [Full summary and execution history](summary.json).
The summary combines deterministic ED, DMRG, thermal-ED and free-DQMC arms with
corrected interacting-DQMC arms. It excludes all original-kernel stochastic
passes and the runtime pilot. The four precision-limited corrected cases use
the registered 16,384-measurement stage, original seeds and unchanged margins.

These are 100 small finite-system observable/parameter comparisons across
repeated reference families. They do not complete the separate 100-idea research
programme or certify every model, size, continuum limit or sign-problem regime.
The sampling defect found in the original atom calculation and all intermediate
results remain public alongside the correction.
