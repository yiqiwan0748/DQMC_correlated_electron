using AltermagnetDQMCTrial
using LinearAlgebra: BLAS

BLAS.set_num_threads(1)

const OUTPUT_DIRECTORY = normpath(joinpath(
    @__DIR__,
    "..",
    "results",
    "dqmc_smoke",
))
const SUMMARY_PATH = joinpath(OUTPUT_DIRECTORY, "validation_summary.csv")
const METADATA_PATH = joinpath(OUTPUT_DIRECTORY, "metadata.txt")

function validate_smoke_summary(summary)
    0 <= summary.acceptance_rate <= 1 ||
        error("invalid acceptance rate for $(summary.case)")
    -1 <= summary.average_sign <= 1 ||
        error("invalid average sign for $(summary.case)")
    0 <= summary.total_density <= 4 ||
        error("invalid total density for $(summary.case)")
    0 <= summary.orbital_x_density <= 2 ||
        error("invalid d_xz density for $(summary.case)")
    0 <= summary.orbital_y_density <= 2 ||
        error("invalid d_yz density for $(summary.case)")
    isapprox(
        summary.total_density,
        summary.orbital_x_density + summary.orbital_y_density;
        atol=2.0e-12,
        rtol=0.0,
    ) || error("orbital densities do not sum to total density for $(summary.case)")
    0 <= summary.double_occupancy_x <= 1 ||
        error("invalid d_xz double occupancy for $(summary.case)")
    0 <= summary.double_occupancy_y <= 1 ||
        error("invalid d_yz double occupancy for $(summary.case)")
    summary.maximum_greens_recompute_error >= 0 ||
        error("negative Green recomputation error for $(summary.case)")
    summary.runtime_seconds >= 0 ||
        error("negative runtime for $(summary.case)")
    numerical_values = (
        summary.acceptance_rate,
        summary.average_sign,
        summary.total_density,
        summary.orbital_x_density,
        summary.orbital_y_density,
        summary.double_occupancy_x,
        summary.double_occupancy_y,
        summary.maximum_greens_recompute_error,
        summary.runtime_seconds,
    )
    all(isfinite, numerical_values) ||
        error("non-finite smoke result for $(summary.case)")
    return nothing
end

cases = (
    (
        name="u0_integration",
        parameters=DQMCParameters(
            L=2,
            beta=1.0,
            dtau=0.1,
            U=0.0,
            mu=0.0,
            stabilization_period=5,
            burnin_sweeps=20,
            measurement_sweeps=50,
            sweeps_between_measurements=1,
            seed=5101,
        ),
    ),
    (
        name="interacting_mu0",
        parameters=DQMCParameters(
            L=2,
            beta=1.0,
            dtau=0.1,
            U=2.0,
            mu=0.0,
            stabilization_period=5,
            burnin_sweeps=20,
            measurement_sweeps=50,
            sweeps_between_measurements=1,
            seed=5102,
        ),
    ),
    (
        name="interacting_first_vhs",
        parameters=DQMCParameters(
            L=2,
            beta=1.0,
            dtau=0.1,
            U=2.0,
            mu=-1.897,
            stabilization_period=5,
            burnin_sweeps=20,
            measurement_sweeps=50,
            sweeps_between_measurements=1,
            seed=5103,
        ),
    ),
)

summaries = map(cases) do case
    summary = run_dqmc_smoke(case.parameters; case_name=case.name)
    validate_smoke_summary(summary)
    return summary
end

columns = (
    :case,
    :L,
    :beta,
    :dtau,
    :Ltau,
    :U,
    :mu,
    :matrix_dimension,
    :burnin_sweeps,
    :measurement_sweeps,
    :sweeps_between_measurements,
    :seed,
    :acceptance_rate,
    :average_sign,
    :total_density,
    :orbital_x_density,
    :orbital_y_density,
    :double_occupancy_x,
    :double_occupancy_y,
    :maximum_greens_recompute_error,
    :runtime_seconds,
)

mkpath(OUTPUT_DIRECTORY)
open(SUMMARY_PATH, "w") do io
    println(io, join(string.(columns), ','))
    for summary in summaries
        println(io, join(
            (string(getproperty(summary, column)) for column in columns),
            ',',
        ))
    end
end

metadata = """
Milestone 5 two-orbital DQMC deterministic smoke validation

Purpose
-------
These deliberately tiny runs validate software behavior. They are not a
parameter scan and support no physical conclusion.

Hamiltonian and chemical potential
----------------------------------
H = sum_sigma c_sigma^dagger K c_sigma
    + U sum_(r,a) (n_up,r,a - 1/2)(n_down,r,a - 1/2)
    - mu sum_(sigma,r,a) n_sigma,r,a.
K is the validated dense real-space hopping matrix and K_mu = K - mu*I, so a
positive mu lowers every one-particle energy. The interaction is centered at
one half exactly as in reference/square_hubbard.jl.

Orbital and index conventions
-----------------------------
Orbital 1 is d_xz and orbital 2 is d_yz. For x,y in 1:L,
site = x + L*(y-1) and p = 2*(site-1) + orbital. One spin sector has
Nsp = 2*L^2 rows. The HS field has shape Nsp by Ltau.

Imaginary-time grid
-------------------
Ltau = beta/dtau must be integral within
64*eps(T)*max(1,abs(beta/dtau)), with zero relative tolerance.

Discrete Hubbard-Stratonovich convention
-----------------------------------------
exp[-dtau*U*(n_up-1/2)*(n_down-1/2)]
  = (1/2)*exp(-dtau*U/4) * sum_(s=+1,-1)
      exp[alpha*s*(n_up-n_down)],
alpha = acosh(exp(dtau*U/2)).
The spin-up potential is V_up = -alpha*s/dtau and the spin-down potential is
V_down = +alpha*s/dtau. Therefore exp(-dtau*V_up) = exp(+alpha*s) and
exp(-dtau*V_down) = exp(-alpha*s). The field-independent factor
(1/2)*exp(-dtau*U/4) for every site-orbital and time slice is omitted from
Metropolis ratios because it cancels between configurations.

Propagator and Green conventions
--------------------------------
Dense exact kinetic exponentials are used, checkerboard = false, and the
asymmetric convention is
B_sigma,l = exp(-dtau*V_sigma,l) * exp(-dtau*K_mu).
The full product is B_sigma,Ltau * ... * B_sigma,1. The equal-time convention
is G_sigma[i,j] = <c_sigma,i c_sigma,j^dagger> and
G_sigma(0,0) = inv(I + B_sigma,Ltau*...*B_sigma,1).

Determinants, signs, and updates
-------------------------------
JDQMCFramework stores log(abs(det(G_sigma))) and sign(det(G_sigma)). The
fermion configuration weight is det(G_up^-1)*det(G_down^-1). Updates sample
its absolute value with probability min(1,abs(R_up*R_down)); the observable
sign is sign(det(G_up^-1)*det(G_down^-1)). Reported densities and double
occupancies are sign reweighted. A free positive weight does not imply that
the interacting two-orbital model is sign-problem free.

Stabilization and errors
------------------------
The JDQMCFramework stabilization period is 5 slices for all three runs.
maximum_greens_recompute_error is the largest elementwise correction reported
by periodic stabilization or the final full Green recomputation. Runtime
starts after state construction and covers burn-in, measurement-separation
sweeps, measurements, and the final recomputation. Acceptance includes both
burn-in and measurement-separation sweeps.

Smoke cases and seeds
---------------------
u0_integration:             L=2, beta=1, dtau=0.1, U=0, mu=0, seed=5101
interacting_mu0:            L=2, beta=1, dtau=0.1, U=2, mu=0, seed=5102
interacting_first_vhs:      L=2, beta=1, dtau=0.1, U=2, mu=-1.897, seed=5103
Each case uses 20 burn-in sweeps, 50 measurements, and one full sweep between
measurements. Fixed seeds reproduce all sampled values; wall-clock runtime is
machine- and load-dependent.

Scope
-----
No spin structure factor, altermagnetic susceptibility, transport observable,
production jackknife analysis, or parameter scan is included.
"""

open(METADATA_PATH, "w") do io
    write(io, metadata)
end

for summary in summaries
    println(
        "SMOKE_RESULT case=", summary.case,
        " acceptance_rate=", summary.acceptance_rate,
        " average_sign=", summary.average_sign,
        " total_density=", summary.total_density,
        " runtime_seconds=", summary.runtime_seconds,
        " maximum_greens_recompute_error=",
        summary.maximum_greens_recompute_error,
    )
end
println("Wrote ", SUMMARY_PATH)
println("Wrote ", METADATA_PATH)
