# ALF eligibility v2: preserved setup failure

Job 6458532 failed after 1m02s, before the observer rebuild and every configuration/buffer proof. The strict cleanup guard found `Prog/opt_time_dependent_mod.mod`; upstream cleanup had removed the differently named `opt_time_dependent.mod`. The guard stopped the run as designed.

P21 remains **NOT_RUN**. This is a setup failure, with no external statistical comparison or scientific discrepancy. The original-native single-provider linkage preflight passed, which does not establish observer proof completion.

Seven original text receipts/logs were retrieved and verified against `remote_sha256.txt`. `slurm_status.txt` records the read-only scheduler status; `failure_audit.json` closes the text evidence with hashes. No binaries, libraries, symlinks or source-tree copies were retrieved. Prior versions and frozen source plans are preserved. A separate v3 cleanup proposal is reviewed separately.
