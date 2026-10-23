using Printf
using Random
using LinearAlgebra

include(joinpath(@__DIR__, "..", "archive", "validated_old_pipeline", "run_lambda_grid_step_validated.jl"))

function timer_dict()
    keys = [
        :wrapper_alloc,
        :initialize_unequal,
        :propagate_up,
        :propagate_dn,
        :inner_total,
        :shuffle,
        :det_ratio,
        :copy_G00_full,
        :logX_old,
        :tentative_G00_update,
        :logX_new,
        :scalar_acceptance,
        :accepted_extra_rank1,
        :accepted_equal_time_update,
        :reject_restore,
        :stabilize_up,
        :stabilize_dn,
        :other_or_loop_overhead,
    ]
    return Dict(k => 0.0 for k in keys)
end

@inline function add_time!(timers, key, t0)
    timers[key] += 1.0e-9 * (time_ns() - t0)
    return nothing
end

function _local_update_replica_fused_timed!(
    Gup, Gup_00, Gup_τ0, Gup_0τ, logdetGup, sgndetGup, Bup,
    Gdn, Gdn_00, Gdn_τ0, Gdn_0τ, logdetGdn, sgndetGdn, Bdn,
    G2up_00, G2dn_00,
    s, μ, α, Δτ, rng, perm, u, v,
    C, λ,
    timers, counts
)
    t = time_ns()
    shuffle!(rng, perm)
    add_time!(timers, :shuffle, t)

    accepted = 0
    N = length(perm)

    t = time_ns()
    logX_current = logabsX_from_G00(Gup_00, Gdn_00, G2up_00, G2dn_00, C)
    add_time!(timers, :logX_old, t)

    for i in perm
        counts[:proposals] += 1

        t = time_ns()
        ΔVup = -α / Δτ * (-2 * s[i])
        ΔVdn = +α / Δτ * (-2 * s[i])
        Rup, Δup = jdqmcf.local_update_det_ratio(Gup, ΔVup, i, Δτ)
        Rdn, Δdn = jdqmcf.local_update_det_ratio(Gdn, ΔVdn, i, Δτ)
        add_time!(timers, :det_ratio, t)

        t = time_ns()
        Gup_00_old = copy(Gup_00)
        Gdn_00_old = copy(Gdn_00)
        add_time!(timers, :copy_G00_full, t)

        t = time_ns()
        logX_old = logX_current
        add_time!(timers, :logX_old, t)

        t = time_ns()
        local_update_G00!(Gup_00, Gup_τ0, Gup_0τ, Rup, Δup, i, u, v)
        local_update_G00!(Gdn_00, Gdn_τ0, Gdn_0τ, Rdn, Δdn, i, u, v)
        add_time!(timers, :tentative_G00_update, t)

        t = time_ns()
        logX_new = logabsX_from_G00(Gup_00, Gdn_00, G2up_00, G2dn_00, C)
        add_time!(timers, :logX_new, t)

        t = time_ns()
        Rλ = fused_ratio_from_logs(logX_new, logX_old, λ)
        P = min(1.0, abs(Rup * Rdn) * Rλ)
        accept = rand(rng) < P
        add_time!(timers, :scalar_acceptance, t)

        if accept
            accepted += 1
            counts[:accepted] += 1
            s[i] = -s[i]

            t = time_ns()
            local_update_Gτ0!(Gup_τ0, Gup, Rup, Δup, i, u, v)
            local_update_G0τ!(Gup_0τ, Gup, Rup, Δup, i, u, v)
            local_update_Gτ0!(Gdn_τ0, Gdn, Rdn, Δdn, i, u, v)
            local_update_G0τ!(Gdn_0τ, Gdn, Rdn, Δdn, i, u, v)
            add_time!(timers, :accepted_extra_rank1, t)

            t = time_ns()
            logdetGup, sgndetGup = jdqmcf.local_update_greens!(
                Gup, logdetGup, sgndetGup, Bup, Rup, Δup, i, u, v
            )
            logdetGdn, sgndetGdn = jdqmcf.local_update_greens!(
                Gdn, logdetGdn, sgndetGdn, Bdn, Rdn, Δdn, i, u, v
            )
            add_time!(timers, :accepted_equal_time_update, t)

            logX_current = logX_new
        else
            counts[:rejected] += 1

            t = time_ns()
            copyto!(Gup_00, Gup_00_old)
            copyto!(Gdn_00, Gdn_00_old)
            add_time!(timers, :reject_restore, t)
        end
    end

    acceptance_rate = accepted / N
    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, acceptance_rate
end

function local_update_replica_fused_timed!(
    Gup::Matrix{T}, logdetGup, sgndetGup, Bup, fermion_greens_calc_up,
    Gdn::Matrix{T}, logdetGdn, sgndetGdn, Bdn, fermion_greens_calc_dn,
    Gup_00, Gup_τ0, Gup_0τ,
    Gdn_00, Gdn_τ0, Gdn_0τ,
    G2up_00, G2dn_00,
    s, μ, α, Δτ, δG, rng,
    C, λ,
    timers, counts
) where {T<:Number}

    Lτ = length(Bup)
    N = size(Gup, 1)

    t = time_ns()
    u = zeros(T, N)
    v = zeros(T, N)
    perm = collect(1:N)
    add_time!(timers, :wrapper_alloc, t)

    acceptance_rate = 0.0

    t = time_ns()
    jdqmcf.initialize_unequaltime_greens!(Gup_τ0, Gup_0τ, Gup, Gup_00)
    jdqmcf.initialize_unequaltime_greens!(Gdn_τ0, Gdn_0τ, Gdn, Gdn_00)
    add_time!(timers, :initialize_unequal, t)

    for l in fermion_greens_calc_up
        t = time_ns()
        jdqmcf.propagate_unequaltime_greens!(
            Gup_τ0, Gup_0τ, Gup,
            fermion_greens_calc_up, Bup
        )
        add_time!(timers, :propagate_up, t)

        t = time_ns()
        jdqmcf.propagate_unequaltime_greens!(
            Gdn_τ0, Gdn_0τ, Gdn,
            fermion_greens_calc_dn, Bdn
        )
        add_time!(timers, :propagate_dn, t)

        s_l = @view s[:, l]

        t = time_ns()
        logdetGup, sgndetGup,
        logdetGdn, sgndetGdn,
        acceptance_rate_l = _local_update_replica_fused_timed!(
            Gup, Gup_00, Gup_τ0, Gup_0τ, logdetGup, sgndetGup, Bup[l],
            Gdn, Gdn_00, Gdn_τ0, Gdn_0τ, logdetGdn, sgndetGdn, Bdn[l],
            G2up_00, G2dn_00,
            s_l, μ, α, Δτ, rng, perm, u, v,
            C, λ,
            timers, counts
        )
        add_time!(timers, :inner_total, t)

        acceptance_rate += acceptance_rate_l / Lτ

        t = time_ns()
        logdetGup, sgndetGup, δGup, δθup = jdqmcf.stabilize_unequaltime_greens!(
            Gup_τ0, Gup_0τ, Gup, Gup_00,
            logdetGup, sgndetGup,
            fermion_greens_calc_up, Bup,
            update_B̄ = true
        )
        add_time!(timers, :stabilize_up, t)

        t = time_ns()
        logdetGdn, sgndetGdn, δGdn, δθdn = jdqmcf.stabilize_unequaltime_greens!(
            Gdn_τ0, Gdn_0τ, Gdn, Gdn_00,
            logdetGdn, sgndetGdn,
            fermion_greens_calc_dn, Bdn,
            update_B̄ = true
        )
        add_time!(timers, :stabilize_dn, t)

        δG = max(δG, δGup, δGdn)

        iterate(fermion_greens_calc_dn, fermion_greens_calc_up)
    end

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, δG, acceptance_rate
end

function update_replica_fused_timed!(active, passive, δG, C, λ, timers, counts)
    logdetGup, sgndetGup,
    logdetGdn, sgndetGdn,
    δG, acceptance_rate = local_update_replica_fused_timed!(
        active.Gup, active.logdetGup, active.sgndetGup,
        active.Bup, active.fermion_greens_calc_up,

        active.Gdn, active.logdetGdn, active.sgndetGdn,
        active.Bdn, active.fermion_greens_calc_dn,

        active.Gup_00, active.Gup_τ0, active.Gup_0τ,
        active.Gdn_00, active.Gdn_τ0, active.Gdn_0τ,

        passive.Gup_00, passive.Gdn_00,

        active.s, active.μ, active.α, active.Δτ,
        δG, active.rng,

        C, λ,
        timers, counts
    )

    active = merge(active, (
        logdetGup = logdetGup,
        sgndetGup = sgndetGup,
        logdetGdn = logdetGdn,
        sgndetGdn = sgndetGdn,
    ))

    return active, δG, acceptance_rate
end

function run_microtiming(; updates_to_time = 20)
    L = 8
    t = 1.0
    U = 2.0
    μ = 0.0
    β = 32.0
    Δτ = 0.05
    n_stab = 10
    C = collect(1:4)
    λ = 0.5

    # Warm compile path on disposable replicas.
    warm1 = initializer(t, U, μ, β, L, Δτ, 92001; n_stab = n_stab)
    warm2 = initializer(t, U, μ, β, L, Δτ, 92002; n_stab = n_stab)
    warm_timers = timer_dict()
    warm_counts = Dict(:proposals => 0, :accepted => 0, :rejected => 0)
    warm1, _, _ = update_replica_fused_timed!(warm1, warm2, 0.0, C, λ, warm_timers, warm_counts)

    replica1 = initializer(t, U, μ, β, L, Δτ, 93001; n_stab = n_stab)
    replica2 = initializer(t, U, μ, β, L, Δτ, 93002; n_stab = n_stab)

    timers = timer_dict()
    counts = Dict(:proposals => 0, :accepted => 0, :rejected => 0)

    δG = 0.0
    t_wall = time()

    for update in 1:updates_to_time
        replica1, δG, _ = update_replica_fused_timed!(replica1, replica2, δG, C, λ, timers, counts)
        replica2, δG, _ = update_replica_fused_timed!(replica2, replica1, δG, C, λ, timers, counts)
    end

    wall = time() - t_wall
    timed_sum = sum(values(timers)) - timers[:other_or_loop_overhead]
    timers[:other_or_loop_overhead] = max(0.0, wall - timed_sum)

    return (; timers, counts, wall, δG, updates_to_time)
end

function main()
    mkpath("results/processed")

    result = run_microtiming(updates_to_time = 20)
    timers = result.timers
    counts = result.counts
    wall = result.wall

    out = "results/processed/microtime_fused_update_L8_beta32_k2_logXcache.csv"

    ordered = [
        :wrapper_alloc,
        :initialize_unequal,
        :propagate_up,
        :propagate_dn,
        :inner_total,
        :shuffle,
        :det_ratio,
        :copy_G00_full,
        :logX_old,
        :tentative_G00_update,
        :logX_new,
        :scalar_acceptance,
        :accepted_extra_rank1,
        :accepted_equal_time_update,
        :reject_restore,
        :stabilize_up,
        :stabilize_dn,
        :other_or_loop_overhead,
    ]

    open(out, "w") do io
        println(io, "section,seconds,percent_of_wall")
        for key in ordered
            seconds = timers[key]
            percent = wall > 0 ? 100 * seconds / wall : NaN
            @printf(io, "%s,%.9f,%.6f\n", String(key), seconds, percent)
        end
    end

    println("Saved:")
    println(out)
    println()

    @printf("wall_seconds = %.6f\n", wall)
    println("updates_to_time = ", result.updates_to_time)
    println("replica_updates_timed = ", 2 * result.updates_to_time)
    println("proposals = ", counts[:proposals])
    println("accepted = ", counts[:accepted])
    println("rejected = ", counts[:rejected])
    @printf("acceptance = %.6f\n", counts[:accepted] / counts[:proposals])
    @printf("deltaG = %.3e\n", result.δG)
    println()

    println("Top sections:")
    rows = sort(collect(timers), by = x -> x[2], rev = true)
    for (key, seconds) in rows
        percent = wall > 0 ? 100 * seconds / wall : NaN
        @printf("%-32s %12.6f s  %8.3f %%\n", String(key), seconds, percent)
    end
end

main()