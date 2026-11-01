# Distinct package benchmarks

Seven new scientific questions are executable: P01–P06 andP16. They compare actual user-ED outputs with independent QuSpin outputs, with [fixed conventions and tolerances](PROTOCOL.md). TeNPy IDs remain NOT_RUN. The drivers count an entire observable matrix/spectrum as one question.

From the staged `source/distinct_benchmarks/packages` directory onISAAC:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_corrected"
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export DISTINCT_ATTEMPT="$VALIDATION_ROOT/results/distinct_packages_attempt01"
mkdir "$DISTINCT_ATTEMPT"
"$VALIDATION_ROOT/quspin-venv/bin/python" quspin_distinct.py --output-dir "$DISTINCT_ATTEMPT/oracle"
export DISTINCT_TARGET_OUTPUT="$DISTINCT_ATTEMPT/target.toml"
"$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no user_ed_distinct.jl
python3 compare.py --oracle-dir "$DISTINCT_ATTEMPT/oracle" --target "$DISTINCT_TARGET_OUTPUT" --output "$DISTINCT_ATTEMPT/results.json"
```

Julia's target uses standard libraries and the recovered user's ED source. QuSpin uses the existing isolated environment. Results are established only by actual execution receipts. A FAIL remains visible and may reveal a source/adapter problem; no result is asserted here.
