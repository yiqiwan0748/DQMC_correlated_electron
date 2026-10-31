##
using Random
using LinearAlgebra
using Statistics
using Printf
import JDQMCFramework as jdqmcf
import JDQMCMeasurements as jdqmcm
import LatticeUtilities as lu

# =============================================================================
# Parameters used for simulation
# =====================================================================

t=1.0 # hopping amplitude
U=2.0 # on-site interaction
μ=0 # chemical potential
β=4.0 # inverse temperature
L=4 # number of spatial sites for a PBC Hubbard chain
Δτ=0.05 # imaginary time step for Trotter slice
Lτ=jdqmcf.eval_length_imaginary_axis(β,Δτ) # number of imaginary time steps
N_burnin=1000 # number of quilibrium sweeps
N_measurements=10000 # number of measurement sweeps
n_stab=10 # number of sweeps between every stabilization recalculation-from-scratch
N_bins=50 # number of bins
N_binsize=N_measurements/N_bins # number of sweeps for every bin
seed = abs(rand(Int)) # seed used for random number generator
rng = Xoshiro(seed) # random number generator
N_sweeps=1 # number of sweeps between every measurement
λ = 1.0

symmetric = false
checkerboard = false




# ===================================================================
# Lattice configurations
# ===================================================================

unit_cell = lu.UnitCell(lattice_vecs = [[1.]], basis_vecs = [[0.]])

# define size of lattice
lattice = lu.Lattice(L = [L], periodic = [true])

# define bonds/hoppings in x direction
bond_x = lu.Bond(orbitals = (1,1), displacement = [1])

# construct neighbor table
neighbor_table = lu.build_neighbor_table([bond_x], unit_cell, lattice)

# calculate number of sites in lattice
N = lu.nsites(unit_cell, lattice)

# calculate number of bonds in lattice
N_bonds = size(neighbor_table, 2)

(N, N_bonds)

Δτ′ = symmetric ? Δτ/2 : Δτ

if checkerboard
    expnΔτ′K = cb.CheckerboardMatrix(neighbor_table, fill(t, N_bonds), Δτ′)
else
    K = zeros(typeof(t), N, N)
    for bond in 1:N_bonds
        i, j = neighbor_table[1, bond], neighbor_table[2, bond]
        K[i,j] = -t
        K[j,i] = -conj(t)
    end
    expnΔτ′K = exp(-Δτ′*K)
    exppΔτ′K = exp(+Δτ′*K)
end;



# ==================================================================
# Initialization (Green's functions, propogators, etc.) 
# ==================================================================

α = acosh(exp(Δτ*U/2))
s = rand(rng, -1:2:1, N, Lτ)
T_expnΔτK = eltype(t)
T_expnΔτV = typeof(α)

if checkerboard && symmetric

    ## Propagator defined as B[σ,l] = exp{-Δτ⋅K/2}⋅exp{-Δτ⋅V[σ,l]}⋅exp{-Δτ⋅K/2},
    ## where the dense matrix exp{-Δτ⋅K/2} is approximated by the sparse checkerboard matrix.
    Bup = jdqmcf.SymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]
    Bdn = jdqmcf.SymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]

elseif checkerboard && !symmetric

    ## Propagator defined as B[σ,l] = exp{-Δτ⋅V[σ,l]}⋅exp{-Δτ⋅K},
    ## where the dense matrix exp{-Δτ⋅K} is approximated by the sparse checkerboard matrix.
    Bup = jdqmcf.AsymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]
    Bdn = jdqmcf.AsymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]

elseif !checkerboard && symmetric

    ## Propagator defined as B[σ,l] = exp{-Δτ⋅K/2}⋅exp{-Δτ⋅V[σ,l]}⋅exp{-Δτ⋅K/2},
    ## where the dense matrix exp{-Δτ⋅K/2} is exactly calculated.
    Bup = jdqmcf.SymExactPropagator{T_expnΔτK, T_expnΔτV}[]
    Bdn = jdqmcf.SymExactPropagator{T_expnΔτK, T_expnΔτV}[]

elseif !checkerboard && !symmetric

    ## Propagator defined as B[σ,l] = exp{-Δτ⋅V[σ,l]}⋅exp{-Δτ⋅K},
    ## where the dense matrix exp{-Δτ⋅K} is exactly calculated.
    Bup = jdqmcf.AsymExactPropagator{T_expnΔτK, T_expnΔτV}[]
    Bdn = jdqmcf.AsymExactPropagator{T_expnΔτK, T_expnΔτV}[]
end;

for l in 1:Lτ

    ## Get the HS fields associated with the current time-slice l.
    s_l = @view s[:,l]

    ## Calculate the spin-up diagonal exponentiated potential energy
    ## matrix exp{-Δτ⋅V[↑,l]} = exp{-Δτ⋅(-α/Δτ⋅s[i,l]-μ)} = exp{+α⋅s[i,l] + Δτ⋅μ}.
    expnΔτVup = zeros(T_expnΔτV, N)
    @. expnΔτVup = exp(+α * s_l + Δτ*μ)

    ## Calculate the spin-down diagonal exponentiated potential energy
    ## matrix exp{-Δτ⋅V[↓,l]} = exp{-Δτ⋅(+α/Δτ⋅s[i,l]-μ)} = exp{-α⋅s[i,l] + Δτ⋅μ}.
    expnΔτVdn = zeros(T_expnΔτV, N)
    @. expnΔτVdn = exp(-α * s_l + Δτ*μ)
    
    ## Initialize spin-up and spin-down propagator matrix for the current time-slice l.
    if checkerboard && symmetric

        push!(Bup, jdqmcf.SymChkbrdPropagator(expnΔτVup, expnΔτ′K))
        push!(Bdn, jdqmcf.SymChkbrdPropagator(expnΔτVdn, expnΔτ′K))

    elseif checkerboard && !symmetric

        push!(Bup, jdqmcf.AsymChkbrdPropagator(expnΔτVup, expnΔτ′K))
        push!(Bdn, jdqmcf.AsymChkbrdPropagator(expnΔτVdn, expnΔτ′K))

    elseif !checkerboard && symmetric

        push!(Bup, jdqmcf.SymExactPropagator(expnΔτVup, expnΔτ′K, exppΔτ′K))
        push!(Bdn, jdqmcf.SymExactPropagator(expnΔτVdn, expnΔτ′K, exppΔτ′K))

    elseif !checkerboard && !symmetric

        push!(Bup, jdqmcf.AsymExactPropagator(expnΔτVup, expnΔτ′K, exppΔτ′K))
        push!(Bdn, jdqmcf.AsymExactPropagator(expnΔτVdn, expnΔτ′K, exppΔτ′K))
    end
end

fermion_greens_calc_up = jdqmcf.FermionGreensCalculator(Bup, β, Δτ, n_stab)
fermion_greens_calc_dn = jdqmcf.FermionGreensCalculator(Bdn, β, Δτ, n_stab);

Gup = zeros(typeof(t), N, N)
logdetGup, sgndetGup = jdqmcf.calculate_equaltime_greens!(Gup, fermion_greens_calc_up)

Gdn = zeros(typeof(t), N, N)
logdetGdn, sgndetGdn = jdqmcf.calculate_equaltime_greens!(Gdn, fermion_greens_calc_dn);

Gup_τ0 = zero(Gup) # Gup(τ,0)
Gup_0τ = zero(Gup) # Gup(0,τ)
Gup_00 = zero(Gup) # Gup(0,0)
Gdn_τ0 = zero(Gdn) # Gdn(τ,0)
Gdn_0τ = zero(Gdn) # Gdn(0,τ)
Gdn_00 = zero(Gdn); # Gdn(0,0)

avg_sign = zeros(eltype(Gup), N_bins)
avg_X = zeros(eltype(Gup), N_bins)

# ====================================================================
# helper functions
# =====================================================================

function local_update_G00!(
    G00::AbstractMatrix{T},
    Gτ0::AbstractMatrix{T},
    G0τ::AbstractMatrix{T},
    R,
    Δ,
    i::Int,
    u::AbstractVector{T},
    v::AbstractVector{T}
)::Nothing where {T<:Number}

    copyto!(u, @view G0τ[:, i])

    # vᵀ = G(τ,0)[i, :]
    vt = adjoint(v)
    copyto!(vt, @view Gτ0[i, :])

    # G(0,0) ← G(0,0) + (Δ/R) G(0,τ)[:,i] G(τ,0)[i,:]
    BLAS.ger!(Δ / R, u, v, G00)

    return nothing
end

function local_update_Gτ0!(
    Gτ0::AbstractMatrix{T},
    Gττ::AbstractMatrix{T},
    R,
    Δ,
    i::Int,
    u::AbstractVector{T},
    v::AbstractVector{T}
)::Nothing where {T<:Number}

    copyto!(u, @view Gττ[:, i])

    # vᵀ = G(τ,0)[i, :]
    vt = adjoint(v)
    copyto!(vt, @view Gτ0[i, :])

    BLAS.ger!(Δ / R, u, v, Gτ0)

    return nothing
end

function local_update_G0τ!(
    G0τ::AbstractMatrix{T},
    Gττ::AbstractMatrix{T},
    R,
    Δ,
    i::Int,
    u::AbstractVector{T},
    v::AbstractVector{T}
)::Nothing where {T<:Number}

    copyto!(u, @view G0τ[:, i])

    # vᵀ = G(τ,τ)[i,:] - I[i,:]
    vt = adjoint(v)
    copyto!(vt, @view Gττ[i, :])
    v[i] -= one(T)

    BLAS.ger!(Δ / R, u, v, G0τ)

    return nothing
end

function grover_matrix(G1::AbstractMatrix, G2::AbstractMatrix, C::Vector{Int})
    T = promote_type(eltype(G1), eltype(G2))
    nC = length(C)

    IC = Matrix{T}(I, nC, nC)

    if nC == 0
        return IC
    end

    G1C = Matrix{T}(G1[C, C])
    G2C = Matrix{T}(G2[C, C])

    gC = G1C * G2C + (IC - G1C) * (IC - G2C)

    return gC
end
function logabsX_from_G00(Gup_00, Gdn_00, G2up_00, G2dn_00, C)
    gCup = grover_matrix(Gup_00, G2up_00, C)
    gCdn = grover_matrix(Gdn_00, G2dn_00, C)

    X = det(gCup) * det(gCdn)
    aX = abs(X)

    return aX == 0 ? -Inf : log(aX)
end

function fused_ratio_from_logs(logX_new, logX_old, λ)
    if λ == 0
        return 1.0
    elseif logX_new == -Inf && logX_old == -Inf
        return 1.0
    elseif logX_new == -Inf
        return 0.0
    elseif logX_old == -Inf
        return Inf
    else
        logR = λ * (logX_new - logX_old)
        return exp(clamp(logR, -700.0, 700.0))
    end
end

function initializer(t, U, μ, β, L, Δτ, seed; n_stab=10)
    Lτ=jdqmcf.eval_length_imaginary_axis(β,Δτ)
    rng = Xoshiro(seed) 
    symmetric = false
    checkerboard = false
    unit_cell = lu.UnitCell(lattice_vecs = [[1.]], basis_vecs = [[0.]])
    lattice = lu.Lattice(L = [L], periodic = [true])

    bond_x = lu.Bond(orbitals = (1,1), displacement = [1])

    neighbor_table = lu.build_neighbor_table([bond_x], unit_cell, lattice)

    N = lu.nsites(unit_cell, lattice)

    N_bonds = size(neighbor_table, 2)

    Δτ′ = symmetric ? Δτ/2 : Δτ 

    if checkerboard
        expnΔτ′K = cb.CheckerboardMatrix(neighbor_table, fill(t, N_bonds), Δτ′)
    else
        K = zeros(typeof(t), N, N)
        for bond in 1:N_bonds
            i, j = neighbor_table[1, bond], neighbor_table[2, bond]
            K[i,j] = -t
            K[j,i] = -conj(t)
        end
        expnΔτ′K = exp(-Δτ′*K)
        exppΔτ′K = exp(+Δτ′*K)
    end;

    α = acosh(exp(Δτ*U/2))
    s = rand(rng, -1:2:1, N, Lτ)
    T_expnΔτK = eltype(t)
    T_expnΔτV = typeof(α)

    if checkerboard && symmetric
        
        Bup = jdqmcf.SymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]
        Bdn = jdqmcf.SymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]

    elseif checkerboard && !symmetric
        
        Bup = jdqmcf.AsymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]
        Bdn = jdqmcf.AsymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]

    elseif !checkerboard && symmetric

        Bup = jdqmcf.SymExactPropagator{T_expnΔτK, T_expnΔτV}[]
        Bdn = jdqmcf.SymExactPropagator{T_expnΔτK, T_expnΔτV}[]

    elseif !checkerboard && !symmetric

        Bup = jdqmcf.AsymExactPropagator{T_expnΔτK, T_expnΔτV}[]
        Bdn = jdqmcf.AsymExactPropagator{T_expnΔτK, T_expnΔτV}[]
    end;

    for l in 1:Lτ

        s_l = @view s[:,l]

        expnΔτVup = zeros(T_expnΔτV, N)
        @. expnΔτVup = exp(+α * s_l + Δτ*μ)

        expnΔτVdn = zeros(T_expnΔτV, N)
        @. expnΔτVdn = exp(-α * s_l + Δτ*μ)
        
        if checkerboard && symmetric

            push!(Bup, jdqmcf.SymChkbrdPropagator(expnΔτVup, expnΔτ′K))
            push!(Bdn, jdqmcf.SymChkbrdPropagator(expnΔτVdn, expnΔτ′K))

        elseif checkerboard && !symmetric

            push!(Bup, jdqmcf.AsymChkbrdPropagator(expnΔτVup, expnΔτ′K))
            push!(Bdn, jdqmcf.AsymChkbrdPropagator(expnΔτVdn, expnΔτ′K))

        elseif !checkerboard && symmetric

            push!(Bup, jdqmcf.SymExactPropagator(expnΔτVup, expnΔτ′K, exppΔτ′K))
            push!(Bdn, jdqmcf.SymExactPropagator(expnΔτVdn, expnΔτ′K, exppΔτ′K))

        elseif !checkerboard && !symmetric

            push!(Bup, jdqmcf.AsymExactPropagator(expnΔτVup, expnΔτ′K, exppΔτ′K))
            push!(Bdn, jdqmcf.AsymExactPropagator(expnΔτVdn, expnΔτ′K, exppΔτ′K))
        end
    end

    fermion_greens_calc_up = jdqmcf.FermionGreensCalculator(Bup, β, Δτ, n_stab)
    fermion_greens_calc_dn = jdqmcf.FermionGreensCalculator(Bdn, β, Δτ, n_stab);

    Gup = zeros(typeof(t), N, N)
    logdetGup, sgndetGup = jdqmcf.calculate_equaltime_greens!(Gup, fermion_greens_calc_up)

    Gdn = zeros(typeof(t), N, N)
    logdetGdn, sgndetGdn = jdqmcf.calculate_equaltime_greens!(Gdn, fermion_greens_calc_dn);

    Gup_τ0 = zero(Gup) 
    Gup_0τ = zero(Gup) 
    Gup_00 = copy(Gup) 
    Gdn_τ0 = zero(Gdn) 
    Gdn_0τ = zero(Gdn) 
    Gdn_00 = copy(Gdn);

    jdqmcf.initialize_unequaltime_greens!(Gup_τ0,Gup_0τ,Gup,Gup_00)
    jdqmcf.initialize_unequaltime_greens!(Gdn_τ0,Gdn_0τ,Gdn,Gdn_00)

return (
    t = t,
    U = U,
    μ = μ,
    β = β,
    L = L,
    Δτ = Δτ,
    Lτ = Lτ,
    n_stab = n_stab,
    α=α,

    rng = rng,
    s = s,

    unit_cell = unit_cell,
    lattice = lattice,
    neighbor_table = neighbor_table,

    Bup = Bup,
    Bdn = Bdn,
    fermion_greens_calc_up = fermion_greens_calc_up,
    fermion_greens_calc_dn = fermion_greens_calc_dn,

    Gup = Gup,
    logdetGup = logdetGup,
    sgndetGup = sgndetGup,
    Gup_τ0 = Gup_τ0,
    Gup_0τ = Gup_0τ,
    Gup_00 = Gup_00,

    Gdn = Gdn,
    logdetGdn = logdetGdn,
    sgndetGdn = sgndetGdn,
    Gdn_τ0 = Gdn_τ0,
    Gdn_0τ = Gdn_0τ,
    Gdn_00 = Gdn_00
)

end

function _local_update_replica!(
    Gup, Gup_00, Gup_τ0, Gup_0τ, logdetGup, sgndetGup, Bup,
    Gdn, Gdn_00, Gdn_τ0, Gdn_0τ, logdetGdn, sgndetGdn, Bdn,
    s, μ, α, Δτ, rng, perm, u, v
)

    shuffle!(rng, perm)

    accepted = 0
    N = length(perm)

    for i in perm
        ΔVup = -α / Δτ * (-2 * s[i])
        ΔVdn = +α / Δτ * (-2 * s[i])

        Rup, Δup = jdqmcf.local_update_det_ratio(Gup, ΔVup, i, Δτ)
        Rdn, Δdn = jdqmcf.local_update_det_ratio(Gdn, ΔVdn, i, Δτ)

        P = min(1.0, abs(Rup * Rdn))

        if rand(rng) < P
            accepted += 1
            s[i] = -s[i]

            local_update_G00!(Gup_00, Gup_τ0, Gup_0τ, Rup, Δup, i, u, v)
            local_update_Gτ0!(Gup_τ0, Gup, Rup, Δup, i, u, v)
            local_update_G0τ!(Gup_0τ, Gup, Rup, Δup, i, u, v)

            logdetGup, sgndetGup = jdqmcf.local_update_greens!(
                Gup, logdetGup, sgndetGup, Bup, Rup, Δup, i, u, v
            )

            local_update_G00!(Gdn_00, Gdn_τ0, Gdn_0τ, Rdn, Δdn, i, u, v)
            local_update_Gτ0!(Gdn_τ0, Gdn, Rdn, Δdn, i, u, v)
            local_update_G0τ!(Gdn_0τ, Gdn, Rdn, Δdn, i, u, v)

            logdetGdn, sgndetGdn = jdqmcf.local_update_greens!(
                Gdn, logdetGdn, sgndetGdn, Bdn, Rdn, Δdn, i, u, v
            )
        end
    end

    acceptance_rate = accepted / N

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, acceptance_rate
end

function _local_update_replica_fused!(
    Gup, Gup_00, Gup_τ0, Gup_0τ, logdetGup, sgndetGup, Bup,
    Gdn, Gdn_00, Gdn_τ0, Gdn_0τ, logdetGdn, sgndetGdn, Bdn,
    G2up_00, G2dn_00,
    s, μ, α, Δτ, rng, perm, u, v,
    C, λ
)

    shuffle!(rng, perm)

    accepted = 0
    N = length(perm)

    # Cache the current Grover log weight for this time slice.
    # It changes only when the tentative G00 update is accepted.
    logX_current = logabsX_from_G00(Gup_00, Gdn_00, G2up_00, G2dn_00, C)

    for i in perm
        ΔVup = -α / Δτ * (-2 * s[i])
        ΔVdn = +α / Δτ * (-2 * s[i])

        Rup, Δup = jdqmcf.local_update_det_ratio(Gup, ΔVup, i, Δτ)
        Rdn, Δdn = jdqmcf.local_update_det_ratio(Gdn, ΔVdn, i, Δτ)

        Gup_00_old = copy(Gup_00)
        Gdn_00_old = copy(Gdn_00)

        logX_old = logX_current

        # Tentatively update G00 only.
        local_update_G00!(Gup_00, Gup_τ0, Gup_0τ, Rup, Δup, i, u, v)
        local_update_G00!(Gdn_00, Gdn_τ0, Gdn_0τ, Rdn, Δdn, i, u, v)

        logX_new = logabsX_from_G00(Gup_00, Gdn_00, G2up_00, G2dn_00, C)

        Rλ = fused_ratio_from_logs(logX_new, logX_old, λ)

        P = min(1.0, abs(Rup * Rdn) * Rλ)

        if rand(rng) < P
            accepted += 1
            s[i] = -s[i]

            local_update_Gτ0!(Gup_τ0, Gup, Rup, Δup, i, u, v)
            local_update_G0τ!(Gup_0τ, Gup, Rup, Δup, i, u, v)

            logdetGup, sgndetGup = jdqmcf.local_update_greens!(
                Gup, logdetGup, sgndetGup, Bup, Rup, Δup, i, u, v
            )

            local_update_Gτ0!(Gdn_τ0, Gdn, Rdn, Δdn, i, u, v)
            local_update_G0τ!(Gdn_0τ, Gdn, Rdn, Δdn, i, u, v)

            logdetGdn, sgndetGdn = jdqmcf.local_update_greens!(
                Gdn, logdetGdn, sgndetGdn, Bdn, Rdn, Δdn, i, u, v
            )

            logX_current = logX_new
        else
            copyto!(Gup_00, Gup_00_old)
            copyto!(Gdn_00, Gdn_00_old)
        end
    end

    acceptance_rate = accepted / N

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, acceptance_rate
end



function local_update_replica!(
    Gup::Matrix{T}, logdetGup, sgndetGup, Bup, fermion_greens_calc_up,
    Gdn::Matrix{T}, logdetGdn, sgndetGdn, Bdn, fermion_greens_calc_dn,
    Gup_00, Gup_τ0, Gup_0τ,
    Gdn_00, Gdn_τ0, Gdn_0τ,
    s, μ, α, Δτ, δG, rng
) where {T<:Number}

    Lτ = length(Bup)
    N = size(Gup, 1)

    u = zeros(T, N)
    v = zeros(T, N)
    perm = collect(1:N)

    acceptance_rate = 0.0

    jdqmcf.initialize_unequaltime_greens!(Gup_τ0, Gup_0τ, Gup, Gup_00)
    jdqmcf.initialize_unequaltime_greens!(Gdn_τ0, Gdn_0τ, Gdn, Gdn_00)

    for l in fermion_greens_calc_up
        jdqmcf.propagate_unequaltime_greens!(
            Gup_τ0, Gup_0τ, Gup,
            fermion_greens_calc_up, Bup
        )

        jdqmcf.propagate_unequaltime_greens!(
            Gdn_τ0, Gdn_0τ, Gdn,
            fermion_greens_calc_dn, Bdn
        )

        s_l = @view s[:, l]

        logdetGup, sgndetGup,
        logdetGdn, sgndetGdn,
        acceptance_rate_l = _local_update_replica!(
            Gup, Gup_00, Gup_τ0, Gup_0τ, logdetGup, sgndetGup, Bup[l],
            Gdn, Gdn_00, Gdn_τ0, Gdn_0τ, logdetGdn, sgndetGdn, Bdn[l],
            s_l, μ, α, Δτ, rng, perm, u, v
        )

        acceptance_rate += acceptance_rate_l / Lτ

        logdetGup, sgndetGup, δGup, δθup = jdqmcf.stabilize_unequaltime_greens!(
            Gup_τ0, Gup_0τ, Gup, Gup_00,
            logdetGup, sgndetGup,
            fermion_greens_calc_up, Bup,
            update_B̄ = true
        )

        logdetGdn, sgndetGdn, δGdn, δθdn = jdqmcf.stabilize_unequaltime_greens!(
            Gdn_τ0, Gdn_0τ, Gdn, Gdn_00,
            logdetGdn, sgndetGdn,
            fermion_greens_calc_dn, Bdn,
            update_B̄ = true
        )

        δG = max(δG, δGup, δGdn)

        iterate(fermion_greens_calc_dn, fermion_greens_calc_up)
    end

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, δG, acceptance_rate
end

function local_update_replica_fused!(
    Gup::Matrix{T}, logdetGup, sgndetGup, Bup, fermion_greens_calc_up,
    Gdn::Matrix{T}, logdetGdn, sgndetGdn, Bdn, fermion_greens_calc_dn,
    Gup_00, Gup_τ0, Gup_0τ,
    Gdn_00, Gdn_τ0, Gdn_0τ,
    G2up_00, G2dn_00,
    s, μ, α, Δτ, δG, rng,
    C, λ
) where {T<:Number}

    Lτ = length(Bup)
    N = size(Gup, 1)

    u = zeros(T, N)
    v = zeros(T, N)
    perm = collect(1:N)

    acceptance_rate = 0.0

    jdqmcf.initialize_unequaltime_greens!(Gup_τ0, Gup_0τ, Gup, Gup_00)
    jdqmcf.initialize_unequaltime_greens!(Gdn_τ0, Gdn_0τ, Gdn, Gdn_00)

    for l in fermion_greens_calc_up
        jdqmcf.propagate_unequaltime_greens!(
            Gup_τ0, Gup_0τ, Gup,
            fermion_greens_calc_up, Bup
        )

        jdqmcf.propagate_unequaltime_greens!(
            Gdn_τ0, Gdn_0τ, Gdn,
            fermion_greens_calc_dn, Bdn
        )

        s_l = @view s[:, l]

        logdetGup, sgndetGup,
        logdetGdn, sgndetGdn,
        acceptance_rate_l = _local_update_replica_fused!(
            Gup, Gup_00, Gup_τ0, Gup_0τ, logdetGup, sgndetGup, Bup[l],
            Gdn, Gdn_00, Gdn_τ0, Gdn_0τ, logdetGdn, sgndetGdn, Bdn[l],
            G2up_00, G2dn_00,
            s_l, μ, α, Δτ, rng, perm, u, v,
            C, λ
        )

        acceptance_rate += acceptance_rate_l / Lτ

        logdetGup, sgndetGup, δGup, δθup = jdqmcf.stabilize_unequaltime_greens!(
            Gup_τ0, Gup_0τ, Gup, Gup_00,
            logdetGup, sgndetGup,
            fermion_greens_calc_up, Bup,
            update_B̄ = true
        )

        logdetGdn, sgndetGdn, δGdn, δθdn = jdqmcf.stabilize_unequaltime_greens!(
            Gdn_τ0, Gdn_0τ, Gdn, Gdn_00,
            logdetGdn, sgndetGdn,
            fermion_greens_calc_dn, Bdn,
            update_B̄ = true
        )

        δG = max(δG, δGup, δGdn)

        iterate(fermion_greens_calc_dn, fermion_greens_calc_up)
    end

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, δG, acceptance_rate
end


function update_replica!(replica, δG)
    logdetGup, sgndetGup,
    logdetGdn, sgndetGdn,
    δG, acceptance_rate = local_update_replica!(
        replica.Gup, replica.logdetGup, replica.sgndetGup,
        replica.Bup, replica.fermion_greens_calc_up,

        replica.Gdn, replica.logdetGdn, replica.sgndetGdn,
        replica.Bdn, replica.fermion_greens_calc_dn,

        replica.Gup_00, replica.Gup_τ0, replica.Gup_0τ,
        replica.Gdn_00, replica.Gdn_τ0, replica.Gdn_0τ,

        replica.s, replica.μ, replica.α, replica.Δτ,
        δG, replica.rng
    )

    replica = merge(replica, (
        logdetGup = logdetGup,
        sgndetGup = sgndetGup,
        logdetGdn = logdetGdn,
        sgndetGdn = sgndetGdn,
    ))

    return replica, δG, acceptance_rate
end

function update_replica_fused!(active, passive, δG, C, λ)

    # Exact fallback to the already-validated independent estimator.
    if λ == 0
        return update_replica!(active, δG)
    end

    logdetGup, sgndetGup,
    logdetGdn, sgndetGdn,
    δG, acceptance_rate = local_update_replica_fused!(
        active.Gup, active.logdetGup, active.sgndetGup,
        active.Bup, active.fermion_greens_calc_up,

        active.Gdn, active.logdetGdn, active.sgndetGdn,
        active.Bdn, active.fermion_greens_calc_dn,

        active.Gup_00, active.Gup_τ0, active.Gup_0τ,
        active.Gdn_00, active.Gdn_τ0, active.Gdn_0τ,

        passive.Gup_00, passive.Gdn_00,

        active.s, active.μ, active.α, active.Δτ,
        δG, active.rng,

        C, λ
    )

    active = merge(active, (
        logdetGup = logdetGup,
        sgndetGup = sgndetGup,
        logdetGdn = logdetGdn,
        sgndetGdn = sgndetGdn,
    ))

    return active, δG, acceptance_rate
end

function make_measurements!(
    sgndetGup,
    sgndetGdn,
    bin,
    avg_phase,
    avg_den,
    gCup,
    gCdn,
    avg_num,
    avg_step,
    λ,
    Δλ
)
    @assert 0.0 <= λ <= 1.0
    @assert Δλ > 0.0
    @assert λ + Δλ <= 1.0 + 10eps(Float64)

    # Fermion-determinant sign of the two-replica configuration.
    sgn =
        sign(sgndetGup * sgndetGdn)

    # Total Grover estimator for the fixed subsystem.
    X_raw =
        det(gCup) * det(gCdn)

    if !isfinite(real(X_raw)) ||
       !isfinite(imag(X_raw))
        error(
            "Encountered nonfinite Grover estimator X = ",
            X_raw
        )
    end

    imaginary_tolerance =
        1e-10 * max(1.0, abs(real(X_raw)))

    if abs(imag(X_raw)) > imaginary_tolerance
        error(
            "Grover estimator has a significant imaginary part: X = ",
            X_raw
        )
    end

    X = real(X_raw)

    # Fractional powers require the positive-X sector.
    if X <= 0.0
        error(
            "Power-increment measurement requires X > 0, but X = ",
            X
        )
    end

    logX = log(X)

    # Undo the X^λ factor in the sampled ensemble.
    reweight =
        exp(-λ * logX)

    # Adjacent-step observable:
    # Z(λ + Δλ) / Z(λ) = <X^Δλ>_λ.
    step_weight =
        exp(Δλ * logX)

    # Common λ-ensemble sign normalization.
    avg_phase[bin] += sgn

    # Endpoint diagnostics requested by the postdoc.
    avg_den[bin] +=
        sgn * reweight

    avg_num[bin] +=
        sgn * X * reweight

    # Adjacent grid ratio numerator.
    avg_step[bin] +=
        sgn * step_weight

    return nothing
end

# ===================================================================
# Simulation run
# ===================================================================

function run_simulation!(
    N_burnin,
    N_bins,
    N_binsize,
    N_sweeps,
    t,
    U,
    μ,
    β,
    L,
    Δτ,
    seed1,
    seed2;
    n_stab = 10,
    C = [1, 2],
    λ = 0.0,
    Δλ = 1.0 - λ
)
    @assert 0.0 <= λ < 1.0
    @assert Δλ > 0.0
    @assert λ + Δλ <= 1.0 + 10 * eps(Float64)

    t_start = time()

    # Largest Green-function error corrected by stabilization.
    δG = 0.0

    # Local-update acceptance accumulator.
    acceptance_rate = 0.0

    replica1 =
        initializer(
            t, U, μ, β, L, Δτ, seed1;
            n_stab = n_stab
        )

    replica2 =
        initializer(
            t, U, μ, β, L, Δτ, seed2;
            n_stab = n_stab
        )

    # One averaged value of each estimator is stored per bin.
    avg_phase = zeros(Float64, N_bins)
    avg_den = zeros(Float64, N_bins)
    avg_num = zeros(Float64, N_bins)
    avg_step = zeros(Float64, N_bins)

    # --------------------------------------------------------
    # Burn-in / thermalization
    # --------------------------------------------------------

    for update in 1:N_burnin
        replica1, δG, acc1 =
            update_replica_fused!(
                replica1,
                replica2,
                δG,
                C,
                λ
            )

        replica2, δG, acc2 =
            update_replica_fused!(
                replica2,
                replica1,
                δG,
                C,
                λ
            )

        acceptance_rate +=
            (acc1 + acc2) / 2
    end

    # --------------------------------------------------------
    # Binned production measurements
    # --------------------------------------------------------

    for bin in 1:N_bins
        for measurement in 1:N_binsize

            for sweep in 1:N_sweeps
                replica1, δG, acc1 =
                    update_replica_fused!(
                        replica1,
                        replica2,
                        δG,
                        C,
                        λ
                    )

                replica2, δG, acc2 =
                    update_replica_fused!(
                        replica2,
                        replica1,
                        δG,
                        C,
                        λ
                    )

                acceptance_rate +=
                    (acc1 + acc2) / 2
            end

            gCup =
                grover_matrix(
                    replica1.Gup_00,
                    replica2.Gup_00,
                    C
                )

            gCdn =
                grover_matrix(
                    replica1.Gdn_00,
                    replica2.Gdn_00,
                    C
                )

            make_measurements!(
                replica1.sgndetGup *
                    replica2.sgndetGup,
                replica1.sgndetGdn *
                    replica2.sgndetGdn,
                bin,
                avg_phase,
                avg_den,
                gCup,
                gCdn,
                avg_num,
                avg_step,
                λ,
                Δλ
            )
        end

        avg_phase[bin] /= N_binsize
        avg_den[bin] /= N_binsize
        avg_num[bin] /= N_binsize
        avg_step[bin] /= N_binsize
    end

    N_update_rounds =
        N_burnin +
        N_bins * N_binsize * N_sweeps

    acceptance_rate /=
        N_update_rounds

    runtime =
        time() - t_start

    return (
        acceptance_rate,
        δG,
        runtime,
        avg_phase,
        avg_den,
        avg_num,
        avg_step
    )
end;

# ===================================================================
# Test run
# ==================================================================

#=seed1 = 1234
seed2 = 5678

acceptance_rate, δG, runtime, avg_sign, avg_X = run_simulation!(
    N_burnin, N_bins, Int(N_binsize), N_sweeps,
    t, U, μ, β, L, Δτ, seed1, seed2;
    n_stab = n_stab,
    C = [1, 2],
    λ = λ
)
println("Acceptance Rate = ", acceptance_rate)
println("Largest Numerical Error = ", δG)
println("Simulation Run Time (sec) = ", runtime)
println("Mean sign = ", mean(avg_sign))
println("Mean X = ", mean(avg_X))
println("S2 estimate = ", -log(abs(mean(avg_X) / mean(avg_sign))))
# ===================================================================
# Data process and visualization
# ===================================================================
# Raw sign average
sign_avg, sign_err = jdqmcm.jackknife(identity, avg_sign)
println("<sign> = ", real(sign_avg), " +/- ", sign_err)

# Raw numerator: <sign * X>
num_avg, num_err = jdqmcm.jackknife(identity, avg_X)
println("<sign * X> = ", real(num_avg), " +/- ", num_err)

# Physical ratio: <sign * X> / <sign>
Xratio_avg, Xratio_err = jdqmcm.jackknife(
    (X, S) -> real(X / S),
    avg_X,
    avg_sign
)
println("Tr rho_A^2 estimator = ", Xratio_avg, " +/- ", Xratio_err)

# Rényi entropy: S2 = -log(<sign * X>/<sign>)
S2_avg, S2_err = jdqmcm.jackknife(
    (X, S) -> -log(abs(X / S)),
    avg_X,
    avg_sign
)
println("S2 = ", real(S2_avg), " +/- ", S2_err)=#

# ===================================================================
# Other stuff
# =====================================================================

# ===================================================================
# Test identity: <X>_0 = 1 / <1/X>_1
# ==================================================================
#=
seed1 = 1234
seed2 = 5678
C_test = [1, 2]

println()
println("========================================")
println("Testing identity: <X>_0 = 1 / <1/X>_1")
println("C = ", C_test)
println("========================================")

# λ = 0 run: directly estimates <X>_0
acc0, dG0, rt0, avg_sign0, avg_X0 = run_simulation!(
    N_burnin, N_bins, Int(N_binsize), N_sweeps,
    t, U, μ, β, L, Δτ, seed1, seed2;
    n_stab = n_stab,
    C = C_test,
    λ = 0.0
)

X0_avg, X0_err = jdqmcm.jackknife(
    (X, S) -> real(X / S),
    avg_X0,
    avg_sign0
)

println()
println("λ = 0 result")
println("Acceptance Rate = ", acc0)
println("δG = ", dG0)
println("<X>_0 = ", X0_avg, " +/- ", X0_err)

# λ = 1 run: denominator estimates <1/X>_1
acc1, dG1, rt1, avg_sign1, avg_X1 = run_simulation!(
    N_burnin, N_bins, Int(N_binsize), N_sweeps,
    t, U, μ, β, L, Δτ, seed1, seed2;
    n_stab = n_stab,
    C = C_test,
    λ = 1.0
)

invX1_avg, invX1_err = jdqmcm.jackknife(identity, avg_sign1)

X_from_lambda1_avg, X_from_lambda1_err = jdqmcm.jackknife(
    S -> real(1 / S),
    avg_sign1
)

println()
println("λ = 1 result")
println("Acceptance Rate = ", acc1)
println("δG = ", dG1)
println("<1/X>_1 = ", invX1_avg, " +/- ", invX1_err)
println("1 / <1/X>_1 = ", X_from_lambda1_avg, " +/- ", X_from_lambda1_err)

println()
println("Identity check")
println("<X>_0              = ", X0_avg, " +/- ", X0_err)
println("1 / <1/X>_1        = ", X_from_lambda1_avg, " +/- ", X_from_lambda1_err)
println("difference          = ", X0_avg - X_from_lambda1_avg)
println("relative difference = ", (X0_avg - X_from_lambda1_avg) / X0_avg)
=#
# ===================================================================

# ===================================================================
# Power-increment ratio processing
# ===================================================================

using Printf

function power_increment_statistics(
    avg_phase,
    avg_den,
    avg_num
)
    phase, phase_err =
        jdqmcm.jackknife(
            identity,
            avg_phase
        )

    den, den_err =
        jdqmcm.jackknife(
            identity,
            avg_den
        )

    num, num_err =
        jdqmcm.jackknife(
            identity,
            avg_num
        )

    # Z₁ / Z_λ
    Z1_over_Zlambda, Z1_over_Zlambda_err =
        jdqmcm.jackknife(
            (num_bin, phase_bin) ->
                real(num_bin / phase_bin),
            avg_num,
            avg_phase
        )

    # Z_λ / Z₀
    Zlambda_over_Z0, Zlambda_over_Z0_err =
        jdqmcm.jackknife(
            (phase_bin, den_bin) ->
                real(phase_bin / den_bin),
            avg_phase,
            avg_den
        )

    # Z₁ / Z₀
    Z1_over_Z0, Z1_over_Z0_err =
        jdqmcm.jackknife(
            (num_bin, den_bin) ->
                real(num_bin / den_bin),
            avg_num,
            avg_den
        )

    # S₂ = -log|Z₁/Z₀|
    S2, S2_err =
        jdqmcm.jackknife(
            (num_bin, den_bin) ->
                -log(abs(num_bin / den_bin)),
            avg_num,
            avg_den
        )

    return (
        phase = real(phase),
        phase_err = phase_err,
        den = real(den),
        den_err = den_err,
        num = real(num),
        num_err = num_err,
        Z1_over_Zlambda = real(Z1_over_Zlambda),
        Z1_over_Zlambda_err = Z1_over_Zlambda_err,
        Zlambda_over_Z0 = real(Zlambda_over_Z0),
        Zlambda_over_Z0_err = Zlambda_over_Z0_err,
        Z1_over_Z0 = real(Z1_over_Z0),
        Z1_over_Z0_err = Z1_over_Z0_err,
        S2 = real(S2),
        S2_err = S2_err
    )
end


# ===================================================================
# CSV summary-row writer
# ===================================================================

function write_power_increment_summary_row!(
    io;
    run_label,
    seed1,
    seed2,
    λ,
    nA,
    S2_ED,
    stats,
    acceptance,
    δG,
    runtime
)
    ΔS2 = stats.S2 - S2_ED

    z =
        stats.S2_err == 0 ?
        NaN :
        ΔS2 / stats.S2_err

    fields = [
        string(run_label),
        string(seed1),
        string(seed2),
        @sprintf("%.3f", λ),
        string(nA),
        @sprintf("%.8f", S2_ED),

        @sprintf("%.8f", stats.Z1_over_Zlambda),
        @sprintf("%.8f", stats.Z1_over_Zlambda_err),

        @sprintf("%.8f", stats.Zlambda_over_Z0),
        @sprintf("%.8f", stats.Zlambda_over_Z0_err),

        @sprintf("%.8f", stats.Z1_over_Z0),
        @sprintf("%.8f", stats.Z1_over_Z0_err),

        @sprintf("%.8f", stats.S2),
        @sprintf("%.8f", stats.S2_err),
        @sprintf("%+.8f", ΔS2),
        @sprintf("%+.4f", z),

        @sprintf("%.8f", stats.phase),
        @sprintf("%.8f", stats.phase_err),

        @sprintf("%.8f", stats.den),
        @sprintf("%.8f", stats.den_err),

        @sprintf("%.8f", stats.num),
        @sprintf("%.8f", stats.num_err),

        isfinite(acceptance) ?
            @sprintf("%.6f", acceptance) :
            "",

        isfinite(δG) ?
            @sprintf("%.3e", δG) :
            "",

        isfinite(runtime) ?
            @sprintf("%.3f", runtime / 3600) :
            ""
    ]

    println(io, join(fields, ","))

    return nothing
end


# ===================================================================
# Run one independent seed pair at a chosen lambda
# ===================================================================


# ===================================================================
# Adjacent lambda-step ratio processing
# ===================================================================

function adjacent_ratio_statistics(
    avg_step,
    avg_phase
)
    step_ratio, step_ratio_err =
        jdqmcm.jackknife(
            (step_bin, phase_bin) ->
                real(step_bin / phase_bin),
            avg_step,
            avg_phase
        )

    entropy_increment, entropy_increment_err =
        jdqmcm.jackknife(
            (step_bin, phase_bin) ->
                -log(abs(step_bin / phase_bin)),
            avg_step,
            avg_phase
        )

    return (
        step_ratio = real(step_ratio),
        step_ratio_err = step_ratio_err,
        entropy_increment =
            real(entropy_increment),
        entropy_increment_err =
            entropy_increment_err
    )
end

function run_power_increment_seed(
    λ,
    seed1,
    seed2;
    Δλ = 1.0 - λ,
    L = 10,
    t = 1.0,
    U = 2.0,
    μ = 0.0,
    β = 40.0,
    Δτ = 0.05,
    n_stab = 10,
    N_sweeps = 1,
    nA = 5,
    N_burnin = 20,
    N_measurements = 40,
    N_bins = 4
)
    @assert 0.0 <= λ < 1.0
    @assert Δλ > 0.0
    @assert λ + Δλ <= 1.0 + 10 * eps(Float64)
    @assert 1 <= nA <= L
    @assert N_burnin >= 0
    @assert N_measurements > 0
    @assert N_bins > 1
    @assert N_measurements % N_bins == 0

    N_binsize =
        div(N_measurements, N_bins)

    C =
        collect(1:nA)

    acceptance,
    δG,
    runtime,
    avg_phase,
    avg_den,
    avg_num,
    avg_step =
        run_simulation!(
            N_burnin,
            N_bins,
            N_binsize,
            N_sweeps,
            t,
            U,
            μ,
            β,
            L,
            Δτ,
            seed1,
            seed2;
            n_stab = n_stab,
            C = C,
            λ = λ,
            Δλ = Δλ
        )

    # Endpoint ratios requested by Jon:
    # Z₁/Z_λ and Z_λ/Z₀.
    stats =
        power_increment_statistics(
            avg_phase,
            avg_den,
            avg_num
        )

    # Adjacent power-increment ratio:
    # Z(λ + Δλ)/Z(λ).
    step_stats =
        adjacent_ratio_statistics(
            avg_step,
            avg_phase
        )

    return (
        λ = λ,
        Δλ = Δλ,
        λ_next = λ + Δλ,
        seed1 = seed1,
        seed2 = seed2,
        C = C,
        acceptance = acceptance,
        δG = δG,
        runtime = runtime,
        avg_phase = avg_phase,
        avg_den = avg_den,
        avg_num = avg_num,
        avg_step = avg_step,
        stats = stats,
        step_stats = step_stats
    )
end


