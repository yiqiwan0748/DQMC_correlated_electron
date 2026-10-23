# X09: certified precision retest

Job 6460656 passed all four actual native DMRG starts at the frozen 16-sweep precision budget, using the original scientific gates. Two manufactured excited-state controls are correctly rejected. The independent stable-residual audit passes 56 checks and the frozen comparator agrees exactly. This adds one X09 pass; X12 already passed and is not counted again.

The original eight-sweep unresolved attempt, its 15 failed detailed gates and cancelled job 6458825 remain in the earlier public record. This is a four-site fixture result, not a claim of generic DMRG convergence. Original physical/native source hashes and all thresholds are unchanged; previously published source freezes are reused.

Reproduce with Python 3.12: `python replay.py --output-dir /tmp/new-x09-analysis`. This checks and extracts all exact raw/member hashes, replays the frozen comparator, then runs the unchanged independent physical-state audit. No solver or sampling is run.
