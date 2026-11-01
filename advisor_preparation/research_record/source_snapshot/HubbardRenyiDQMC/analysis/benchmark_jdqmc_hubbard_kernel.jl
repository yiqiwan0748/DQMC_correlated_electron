using Random
using LinearAlgebra
using Printf

import JDQMCFramework as jdqmcf
import LatticeUtilities as lu

BLAS.set_num_threads(1)

function build_1d_hubbard(; L, t, U, mu, beta, dtau, n_stab, seed)
    rng = Xoshiro(seed)

    unit_cell = lu.UnitCell(
        lattice_vecs = [[1.0]],
        basis_vecs = [[0.0]]
    )

    lattice = lu.Lattice(
        L = [L],
        periodic = [true]
    )

    bond_p = lu.Bond(
        orbitals = (1, 1),
        displacement = [1]
    )

    neighbor_table = lu.build_neighbor_table([bond_p], unit_cell, lattice)

    N = lu.nsites(unit_cell, lattice)
    Ltau = round(Int, beta / dtau)

    K = zeros(Float64, N, N)
    for bond in 1:size(neighbor_table, 2)
        i, j = neighbor_table[1, bond], neighbor_table[2, bond]
        K[i, j] = -t
        K[j, i] = -t
    end

    expn_dtau_K = exp(-dtau * K)
    expp_dtau_K = exp(+dtau * K)

    alpha = acosh(exp(dtau * U / 2))
    s = rand(rng, -1:2:1, N, Ltau)

    Bup = jdqmcf.AsymExactPropagator{Float64, Float64}[]
    Bdn = jdqmcf.AsymExactPropagator{Float64, Float64}[]

    for l in 1:Ltau
        s_l = @view s[:, l]

        expn_dtau_Vup = zeros(Float64, N)
        expn_dtau_Vdn = zeros(Float64, N)

        @. expn_dtau_Vup = exp(+alpha * s_l + dtau * mu)
        @. expn_dtau_Vdn = exp(-alpha * s_l + dtau * mu)

        push!(Bup, jdqmcf.AsymExactPropagator(expn_dtau_Vup, expn_dtau_K, expp_dtau_K))
        push!(Bdn, jdqmcf.AsymExactPropagator(expn_dtau_Vdn, expn_dtau_K, expp_dtau_K))
    end

    calc_up = jdqmcf.FermionGreensCalculator(Bup, beta, dtau, n_stab)
    calc_dn = jdqmcf.FermionGreensCalculator(Bdn, beta, dtau, n_stab)

    Gup = zeros(Float64, N, N)
    Gdn = zeros(Float64, N, N)

    logdetGup, sgndetGup = jdqmcf.calculate_equaltime_greens!(Gup, calc_up)
    logdetGdn, sgndetGdn = jdqmcf.calculate_equaltime_greens!(Gdn, calc_dn)

    return (; rng, N, Ltau, alpha, s, Bup, Bdn, calc_up, calc_dn,
             Gup, Gdn, logdetGup, sgndetGup, logdetGdn, sgndetGdn)
end

function one_timeslice_update!(
    Gup, logdetGup, sgndetGup, Bup_l,
    Gdn, logdetGdn, sgndetGdn, Bdn_l,
    s_l, mu, alpha, dtau, rng, perm, u, v
)
    shuffle!(rng, perm)
    accepted = 0
    N = length(perm)

    for i in perm
        dVup = -alpha / dtau * (-2 * s_l[i])
        dVdn = +alpha / dtau * (-2 * s_l[i])

        Rup, Dup = jdqmcf.local_update_det_ratio(Gup, dVup, i, dtau)
        Rdn, Ddn = jdqmcf.local_update_det_ratio(Gdn, dVdn, i, dtau)

        P = min(1.0, abs(Rup * Rdn))

        if rand(rng) < P
            accepted += 1
            s_l[i] = -s_l[i]

            logdetGup, sgndetGup = jdqmcf.local_update_greens!(
                Gup, logdetGup, sgndetGup, Bup_l, Rup, Dup, i, u, v
            )

            logdetGdn, sgndetGdn = jdqmcf.local_update_greens!(
                Gdn, logdetGdn, sgndetGdn, Bdn_l, Rdn, Ddn, i, u, v
            )
        end
    end

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, accepted / N
end

function local_update_sweep!(state, mu, dtau)
    N = state.N
    A = zeros(Float64, N, N)
    u = zeros(Float64, N)
    v = zeros(Float64, N)
    perm = collect(1:N)

    deltaG = 0.0
    acceptance = 0.0

    logdetGup = state.logdetGup
    sgndetGup = state.sgndetGup
    logdetGdn = state.logdetGdn
    sgndetGdn = state.sgndetGdn

    for l in state.calc_up
        jdqmcf.propagate_equaltime_greens!(state.Gup, state.calc_up, state.Bup)
        jdqmcf.propagate_equaltime_greens!(state.Gdn, state.calc_dn, state.Bdn)

        jdqmcf.partially_wrap_greens_reverse!(state.Gup, state.Bup[l], A)
        jdqmcf.partially_wrap_greens_reverse!(state.Gdn, state.Bdn[l], A)

        s_l = @view state.s[:, l]

        logdetGup, sgndetGup, logdetGdn, sgndetGdn, ac_l = one_timeslice_update!(
            state.Gup, logdetGup, sgndetGup, state.Bup[l],
            state.Gdn, logdetGdn, sgndetGdn, state.Bdn[l],
            s_l, mu, state.alpha, dtau, state.rng, perm, u, v
        )

        acceptance += ac_l / state.Ltau

        jdqmcf.partially_wrap_greens_forward!(state.Gup, state.Bup[l], A)
        jdqmcf.partially_wrap_greens_forward!(state.Gdn, state.Bdn[l], A)

        logdetGup, sgndetGup, dGup, _ = jdqmcf.stabilize_equaltime_greens!(
            state.Gup, logdetGup, sgndetGup, state.calc_up, state.Bup, update_B̄ = true
        )

        logdetGdn, sgndetGdn, dGdn, _ = jdqmcf.stabilize_equaltime_greens!(
            state.Gdn, logdetGdn, sgndetGdn, state.calc_dn, state.Bdn, update_B̄ = true
        )

        deltaG = max(deltaG, dGup, dGdn)

        iterate(state.calc_dn, state.calc_up)
    end

    state = merge(state, (; logdetGup, sgndetGup, logdetGdn, sgndetGdn))
    return state, acceptance, deltaG
end

function run_case(seed)
    L = 8
    t = 1.0
    U = 2.0
    mu = 0.0
    beta = 32.0
    dtau = 0.05
    n_stab = 10
    N_burnin = 500
    N_measurements = 1000
    N_sweeps = 1

    state = build_1d_hubbard(; L, t, U, mu, beta, dtau, n_stab, seed)

    total_sweeps = N_burnin + N_measurements * N_sweeps
    acceptance_sum = 0.0
    deltaG_max = 0.0

    t0 = time()

    for sweep in 1:total_sweeps
        state, ac, deltaG = local_update_sweep!(state, mu, dtau)
        acceptance_sum += ac
        deltaG_max = max(deltaG_max, deltaG)
    end

    wall_seconds = time() - t0
    acceptance = acceptance_sum / total_sweeps
    site_updates = total_sweeps * state.N * state.Ltau
    site_updates_per_second = site_updates / wall_seconds

    return (; seed, L, beta, dtau, n_stab, N_burnin, N_measurements,
            N_sweeps, total_sweeps, site_updates, wall_seconds,
            site_updates_per_second, acceptance, deltaG_max)
end

function main()
    mkpath("results/processed")
    out = "results/processed/speed_compare_jdqmc_kernel_L8_beta32.csv"

    seeds = [90002, 90003, 90004]
    results = map(run_case, seeds)

    open(out, "w") do io
        println(io, "implementation,L,beta,dtau,n_stab,N_burnin,N_measurements,N_sweeps,seed,total_sweeps,site_updates,wall_seconds,site_updates_per_second,acceptance,deltaG_max")
        for r in results
            @printf(io,
                "jdqmc_kernel,%d,%.8f,%.8f,%d,%d,%d,%d,%d,%d,%d,%.6f,%.6f,%.12f,%.12e\n",
                r.L, r.beta, r.dtau, r.n_stab, r.N_burnin, r.N_measurements,
                r.N_sweeps, r.seed, r.total_sweeps, r.site_updates,
                r.wall_seconds, r.site_updates_per_second, r.acceptance, r.deltaG_max
            )
        end
    end

    println("Saved:")
    println(out)

    println()
    println("Results:")
    for r in results
        @printf("seed=%d wall_seconds=%.3f site_updates_per_second=%.3f acceptance=%.6f deltaG_max=%.3e\n",
                r.seed, r.wall_seconds, r.site_updates_per_second, r.acceptance, r.deltaG_max)
    end
end

main()