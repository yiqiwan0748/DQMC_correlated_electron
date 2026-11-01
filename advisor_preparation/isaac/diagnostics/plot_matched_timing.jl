#!/usr/bin/env julia
# Scientific figure only: reads existing timings and never runs a simulation.
# julia --project=/path/to/active/package plot_matched_timing.jl INPUT.toml OUTPUT_PREFIX
ENV["GKSwstype"] = "100"
using TOML, Statistics, SHA, Dates
using Plots

const DIRECT = "without_Grover_cache_direct_fused"
const OPTIMIZED = "with_rank1_Grover_cache"

function finite_positive(x)
    x isa Real && isfinite(x) && x > 0
end

function timing_row(row)
    raw = get(row, "raw_timings", Any[])
    direct = filter(r -> get(r, "variant", "") == DIRECT, raw)
    optimized = filter(r -> get(r, "variant", "") == OPTIMIZED, raw)
    d = [r["node_loop_seconds"] for r in direct if finite_positive(get(r, "node_loop_seconds", NaN))]
    o = [r["node_loop_seconds"] for r in optimized if finite_positive(get(r, "node_loop_seconds", NaN))]
    checks = get(row, "correctness_checks", Any[])
    requested = get(row, "requested_repeats", 0)
    exact_keys = ["sweep_boundary_HS_fields_exact", "sweep_acceptance_rates_exact", "determinant_signs_exact", "next_RNG_draws_exact"]
    complete = get(row, "status", "UNRESOLVED") == "PASS" && requested >= 7 &&
               get(row, "completed_repeats", 0) == requested &&
               length(d) == length(o) == length(checks) == requested &&
               all(c -> get(c, "status", "") == "PASS" && all(k -> get(c, k, false) === true, exact_keys), checks)
    if !complete
        return nothing
    end
    dm, om = median(d), median(o)
    reported = get(row, "ratio_of_median_node_times", NaN)
    isapprox(dm / om, reported; rtol=1e-12, atol=1e-12) || error("Raw median ratio disagrees with receipt at L=$(row["L"])")
    Dict("L" => row["L"], "nA" => get(row, "nA", row["L"] ÷ 2),
         "direct_median_seconds" => dm, "optimized_median_seconds" => om,
         "node_loop_speedup" => dm / om, "paired_repeats" => requested)
end

function main()
    length(ARGS) == 2 || error("Usage: plot_matched_timing.jl INPUT.toml OUTPUT_PREFIX")
    input = abspath(ARGS[1]); prefix = abspath(ARGS[2]); mkpath(dirname(prefix))
    outputs = [prefix * ".png", prefix * ".pdf", prefix * "_plot_receipt.toml"]
    any(isfile, outputs) && error("Choose a new output prefix; existing artifacts are preserved")
    data = TOML.parsefile(input)
    source_rows = get(data, "results", Any[])
    registered = Int.(get(data, "registered_sizes", [8,16,32,64,128]))
    length(unique(registered)) == length(registered) || error("Duplicate registered sizes")
    rows = Any[]; statuses = Dict{String,Any}()
    for L in registered
        matches = filter(r -> get(r, "L", 0) == L, source_rows)
        if isempty(matches)
            statuses[string(L)] = "UNRESOLVED (missing)"
        elseif length(matches) != 1
            statuses[string(L)] = "UNRESOLVED (duplicate records)"
        else
            row = only(matches)
            result = timing_row(row)
            statuses[string(L)] = result === nothing ? get(row, "status", "UNRESOLVED") * " (not a completed production comparison)" : "PASS"
            result === nothing || push!(rows, result)
        end
    end
    sort!(rows; by=r -> r["L"])
    sizes = [r["L"] for r in rows]
    ticks = (registered, string.(registered))
    defaults = (fontfamily="sans-serif", guidefontsize=10, tickfontsize=9, legendfontsize=8,
                framestyle=:box, gridalpha=.18, xscale=:log2, xticks=ticks,
                xlims=(minimum(registered)/1.2, maximum(registered)*1.2),
                xlabel="Chain length L (nA = L/2)")
    left = plot(; defaults..., yscale=:log10, ylabel="Median node-loop elapsed time (s)",
                title="Same workflow, with and without caching", legend=:topleft)
    if !isempty(rows)
        plot!(left, sizes, [r["direct_median_seconds"] for r in rows];
              marker=:circle, markersize=4.5, linewidth=1.8, color="#B34A3C", label="Direct determinant recomputation")
        plot!(left, sizes, [r["optimized_median_seconds"] for r in rows];
              marker=:square, markersize=4.5, linewidth=1.8, color="#245E9B", label="Cached rank-one Grover update")
    else
        plot!(left; ylims=(1e-4,1.0), legend=false)
        annotate!(left, sqrt(minimum(registered)*maximum(registered)), 1e-2,
                  text("No completed passing production sizes",9,:gray40))
    end
    right = plot(; defaults..., ylabel="Direct / optimized node-loop time", legend=false,
                 title="Median speedup (same host and workload)")
    hline!(right, [1.]; color=:gray60, linestyle=:dash, linewidth=1, label="")
    if !isempty(rows)
        ratios = [r["node_loop_speedup"] for r in rows]
        plot!(right, sizes, ratios; marker=:circle, markersize=4.5, linewidth=1.8,
              color="#245E9B", label="")
        low, high = extrema(vcat(1., ratios))
        span = max(high-low, .2)
        plot!(right; ylims=(max(0.,low-.18*span), high+.24*span))
        for r in rows
            annotate!(right, r["L"], r["node_loop_speedup"]+.08*span,
                      text(string(round(r["node_loop_speedup"]; digits=2), "x"),8,:black))
        end
    else
        plot!(right; ylims=(0.,2.))
    end
    status_line = join(["L=$L: $(statuses[string(L)])" for L in registered], " | ")
    # Keep unresolved sizes in the caption and axis ticks; no partial timing can
    # masquerade as a completed comparison. Raw jitter intervals remain in TOML.
    foot = plot(; xlims=(0,1), ylims=(0,1), axis=false, grid=false, legend=false)
    annotate!(foot,.5,.75,text("Two replicas; beta=2, dtau=0.05, U=2, lambda=0.5; seven matched repeats of twelve rounds",9,:black))
    annotate!(foot,.5,.35,text(status_line,8,:gray30))
    annotate!(foot,.5,.02,text("Compilation/initialization excluded; loop GC included. One lambda node, not a complete entropy campaign.",8,:gray30))
    fig = plot(left,right,foot; layout=@layout([a b; c{.17h}]), size=(1200,560),
               dpi=220, margin=5*Plots.mm)
    savefig(fig,outputs[1]); savefig(fig,outputs[2])
    receipt = Dict("created_utc" => string(now(UTC)),
        "input_path" => input, "input_sha256" => bytes2hex(sha256(read(input))),
        "plot_script_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "benchmark_runner_sha256" => get(data,"runner_sha256","MISSING"),
        "registered_sizes" => registered, "retained_statuses" => statuses,
        "plotted_rows" => rows, "png" => outputs[1], "pdf" => outputs[2],
        "julia_version" => string(VERSION), "plots_version" => string(pkgversion(Plots)),
        "scope" => "Completed passing matched node-loop timings only; all registered sizes remain in ticks and caption. This figure adds no simulation data or CV certification.",
        "uncertainty" => "No error bars plotted; paired-bootstrap timing-jitter intervals and all raw repeats remain in source TOML.")
    open(outputs[3],"w") do io; TOML.print(io,receipt); end
    println("Saved ", outputs[1]); println("Saved ", outputs[2])
    println("Completed passing sizes: ", join(sizes,", ")); println(status_line)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
