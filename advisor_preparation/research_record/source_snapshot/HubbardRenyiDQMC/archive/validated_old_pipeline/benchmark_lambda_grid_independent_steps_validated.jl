include("run_lambda_grid_step_validated.jl")

using Statistics
using Printf

# ============================================================
# True power-increment lambda-grid benchmark
#
# Arguments:
#   1: number of increments
#   2: burn-in updates
#   3: measurements per seed and lambda
#   4: bins
#   5: number of independent seed pairs
#
# Example:
#   julia benchmark_lambda_grid.jl 4 1000 4000 20 4
# ============================================================

N_increment =
    length(ARGS) >= 1 ?
    parse(Int, ARGS[1]) :
    4

N_burnin =
    length(ARGS) >= 2 ?
    parse(Int, ARGS[2]) :
    100

N_measurements =
    length(ARGS) >= 3 ?
    parse(Int, ARGS[3]) :
    1000

N_bins =
    length(ARGS) >= 4 ?
    parse(Int, ARGS[4]) :
    20

N_seed_pairs =
    length(ARGS) >= 5 ?
    parse(Int, ARGS[5]) :
    4

@assert N_increment >= 1
@assert N_burnin >= 0
@assert N_measurements > 0
@assert N_bins > 1
@assert N_measurements % N_bins == 0

# ------------------------------------------------------------
# Physical parameters
# ------------------------------------------------------------

L = 10
t = 1.0
U = 2.0
μ = 0.0
β = 40.0
Δτ = 0.05
n_stab = 10
N_sweeps = 1

nA = 5
S2_ED = 1.7165516614333372

available_seed_pairs = [
    (1234, 5678),
    (11234, 15678),
    (21234, 25678),
    (31234, 35678)
]

@assert 1 <= N_seed_pairs <= length(available_seed_pairs)

seed_pairs =
    available_seed_pairs[1:N_seed_pairs]

Δλ =
    1.0 / N_increment

λ_left_values =
    [
        k * Δλ
        for k in 0:(N_increment - 1)
    ]

# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------

run_tag =
    "Ninc$(N_increment)" *
    "_L$(L)_U$(Int(U))_beta$(Int(β))" *
    "_nA$(nA)_nseed$(N_seed_pairs)" *
    "_Nb$(N_burnin)" *
    "_Nm$(N_measurements)" *
    "_Nbin$(N_bins)"

step_csv =
    "dqmc_lambda_grid_$(run_tag)_steps.csv"

bins_csv =
    "dqmc_lambda_grid_$(run_tag)_bins.csv"

total_csv =
    "dqmc_lambda_grid_$(run_tag)_total.csv"

open(step_csv, "w") do io
    println(
        io,
        join(
            [
                "step",
                "run",
                "seed1",
                "seed2",
                "lambda_left",
                "lambda_right",
                "delta_lambda",
                "step_ratio",
                "step_ratio_err",
                "entropy_increment",
                "entropy_increment_err",
                "Z1_over_Zlambda",
                "Z1_over_Zlambda_err",
                "Zlambda_over_Z0",
                "Zlambda_over_Z0_err",
                "acceptance",
                "deltaG",
                "runtime_seconds"
            ],
            ","
        )
    )
end

open(bins_csv, "w") do io
    println(
        io,
        "step,run,bin,lambda_left,lambda_right," *
        "phase,step_estimator,den,num"
    )
end

open(total_csv, "w") do io
    println(
        io,
        "N_increment,delta_lambda,Z1_over_Z0," *
        "Z1_over_Z0_err,S2_DQMC,S2_err,S2_ED," *
        "DeltaS2,z,total_runtime_seconds"
    )
end

# ------------------------------------------------------------
# Run all adjacent lambda intervals
# ------------------------------------------------------------

combined_step_ratios = Float64[]
combined_step_ratio_errors = Float64[]

combined_entropy_increments = Float64[]
combined_entropy_errors = Float64[]

total_runtime = 0.0

println()
println("========================================")
println("Power-increment lambda grid")
println("========================================")
println("N_increment = ", N_increment)
println("delta lambda = ", Δλ)
println("lambda left values = ", λ_left_values)
println("seed pairs = ", length(seed_pairs))
println("========================================")

for (step_id, λ_left) in enumerate(λ_left_values)

    λ_right =
        λ_left + Δλ

    all_phase = Float64[]
    all_step = Float64[]

    step_acceptances = Float64[]
    step_deltaGs = Float64[]
    step_runtimes = Float64[]

    println()
    println("----------------------------------------")
    println(
        "Step ",
        step_id,
        ": ",
        λ_left,
        " -> ",
        λ_right
    )
    println("----------------------------------------")
    flush(stdout)

    for (run_id, (seed1, seed2)) in enumerate(seed_pairs)

        # Use a distinct random stream for every lambda step.
        step_seed1 =
            seed1 + 1_000_000 * (step_id - 1)

        step_seed2 =
            seed2 + 2_000_000 * (step_id - 1)

        result =
            run_power_increment_seed(
                λ_left,
                step_seed1,
                step_seed2;
                Δλ = Δλ,
                L = L,
                t = t,
                U = U,
                μ = μ,
                β = β,
                Δτ = Δτ,
                n_stab = n_stab,
                N_sweeps = N_sweeps,
                nA = nA,
                N_burnin = N_burnin,
                N_measurements = N_measurements,
                N_bins = N_bins
            )

        append!(all_phase, result.avg_phase)
        append!(all_step, result.avg_step)

        push!(step_acceptances, result.acceptance)
        push!(step_deltaGs, result.δG)
        push!(step_runtimes, result.runtime)

        global total_runtime +=
            result.runtime

        step_stats =
            result.step_stats

        endpoint_stats =
            result.stats

        println(
            "run ",
            run_id,
            ": ratio = ",
            step_stats.step_ratio,
            " +/- ",
            step_stats.step_ratio_err,
            ", dS2 = ",
            step_stats.entropy_increment,
            " +/- ",
            step_stats.entropy_increment_err
        )

        open(step_csv, "a") do io
            println(
                io,
                join(
                    [
                        step_id,
                        run_id,
                        step_seed1,
                        step_seed2,
                        @sprintf("%.8f", λ_left),
                        @sprintf("%.8f", λ_right),
                        @sprintf("%.8f", Δλ),
                        @sprintf(
                            "%.12e",
                            step_stats.step_ratio
                        ),
                        @sprintf(
                            "%.12e",
                            step_stats.step_ratio_err
                        ),
                        @sprintf(
                            "%.12e",
                            step_stats.entropy_increment
                        ),
                        @sprintf(
                            "%.12e",
                            step_stats.entropy_increment_err
                        ),
                        @sprintf(
                            "%.12e",
                            endpoint_stats.Z1_over_Zlambda
                        ),
                        @sprintf(
                            "%.12e",
                            endpoint_stats.Z1_over_Zlambda_err
                        ),
                        @sprintf(
                            "%.12e",
                            endpoint_stats.Zlambda_over_Z0
                        ),
                        @sprintf(
                            "%.12e",
                            endpoint_stats.Zlambda_over_Z0_err
                        ),
                        @sprintf(
                            "%.12e",
                            result.acceptance
                        ),
                        @sprintf(
                            "%.12e",
                            result.δG
                        ),
                        @sprintf(
                            "%.6f",
                            result.runtime
                        )
                    ],
                    ","
                )
            )
        end

        open(bins_csv, "a") do io
            for bin in eachindex(result.avg_phase)
                println(
                    io,
                    join(
                        [
                            step_id,
                            run_id,
                            bin,
                            @sprintf("%.8f", λ_left),
                            @sprintf("%.8f", λ_right),
                            @sprintf(
                                "%.12e",
                                result.avg_phase[bin]
                            ),
                            @sprintf(
                                "%.12e",
                                result.avg_step[bin]
                            ),
                            @sprintf(
                                "%.12e",
                                result.avg_den[bin]
                            ),
                            @sprintf(
                                "%.12e",
                                result.avg_num[bin]
                            )
                        ],
                        ","
                    )
                )
            end
        end
    end

    combined_step =
        adjacent_ratio_statistics(
            all_step,
            all_phase
        )

    push!(
        combined_step_ratios,
        combined_step.step_ratio
    )

    push!(
        combined_step_ratio_errors,
        combined_step.step_ratio_err
    )

    push!(
        combined_entropy_increments,
        combined_step.entropy_increment
    )

    push!(
        combined_entropy_errors,
        combined_step.entropy_increment_err
    )

    println()
    println(
        "Combined step ratio = ",
        combined_step.step_ratio,
        " +/- ",
        combined_step.step_ratio_err
    )

    println(
        "Combined entropy increment = ",
        combined_step.entropy_increment,
        " +/- ",
        combined_step.entropy_increment_err
    )

    open(step_csv, "a") do io
        println(
            io,
            join(
                [
                    step_id,
                    "combined",
                    "all",
                    "all",
                    @sprintf("%.8f", λ_left),
                    @sprintf("%.8f", λ_right),
                    @sprintf("%.8f", Δλ),
                    @sprintf(
                        "%.12e",
                        combined_step.step_ratio
                    ),
                    @sprintf(
                        "%.12e",
                        combined_step.step_ratio_err
                    ),
                    @sprintf(
                        "%.12e",
                        combined_step.entropy_increment
                    ),
                    @sprintf(
                        "%.12e",
                        combined_step.entropy_increment_err
                    ),
                    "",
                    "",
                    "",
                    "",
                    @sprintf(
                        "%.12e",
                        mean(step_acceptances)
                    ),
                    @sprintf(
                        "%.12e",
                        maximum(step_deltaGs)
                    ),
                    @sprintf(
                        "%.6f",
                        sum(step_runtimes)
                    )
                ],
                ","
            )
        )
    end
end

# ------------------------------------------------------------
# Multiply adjacent ratios and sum entropy increments
# ------------------------------------------------------------

Z1_over_Z0 =
    prod(combined_step_ratios)

relative_variance =
    sum(
        (
            combined_step_ratio_errors[i] /
            combined_step_ratios[i]
        )^2
        for i in eachindex(combined_step_ratios)
    )

Z1_over_Z0_err =
    abs(Z1_over_Z0) *
    sqrt(relative_variance)

S2_DQMC =
    sum(combined_entropy_increments)

S2_err =
    sqrt(
        sum(
            err^2
            for err in combined_entropy_errors
        )
    )

DeltaS2 =
    S2_DQMC - S2_ED

z =
    S2_err == 0.0 ?
    NaN :
    DeltaS2 / S2_err

println()
println("========================================")
println("Complete lambda-grid result")
println("========================================")

for step_id in eachindex(combined_step_ratios)
    println(
        "step ",
        step_id,
        ": ratio = ",
        combined_step_ratios[step_id],
        " +/- ",
        combined_step_ratio_errors[step_id]
    )
end

println("Z1/Z0 = ", Z1_over_Z0, " +/- ", Z1_over_Z0_err)
println("S2_DQMC = ", S2_DQMC, " +/- ", S2_err)
println("S2_ED = ", S2_ED)
println("DeltaS2 = ", DeltaS2)
println("z = ", z)
println("total runtime seconds = ", total_runtime)

open(total_csv, "a") do io
    println(
        io,
        join(
            [
                N_increment,
                @sprintf("%.8f", Δλ),
                @sprintf("%.12e", Z1_over_Z0),
                @sprintf("%.12e", Z1_over_Z0_err),
                @sprintf("%.12e", S2_DQMC),
                @sprintf("%.12e", S2_err),
                @sprintf("%.12e", S2_ED),
                @sprintf("%.12e", DeltaS2),
                @sprintf("%.12e", z),
                @sprintf("%.6f", total_runtime)
            ],
            ","
        )
    )
end

println()
println("Saved step CSV: ", step_csv)
println("Saved bin CSV: ", bins_csv)
println("Saved total CSV: ", total_csv)


