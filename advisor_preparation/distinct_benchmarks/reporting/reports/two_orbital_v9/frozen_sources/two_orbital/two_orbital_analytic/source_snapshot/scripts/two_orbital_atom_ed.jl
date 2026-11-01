module TwoOrbitalAtomED

using LinearAlgebra
using Serialization

export solve_two_orbital_atom,
       broaden_spectral_peaks,
       save_atom_ed_result

function _annihilation_operator(mode::Int, number_of_modes::Int)
    1 <= mode <= number_of_modes || throw(ArgumentError("invalid mode"))
    dimension = 1 << number_of_modes
    operator = zeros(Float64, dimension, dimension)
    bit = 1 << (mode - 1)
    lower_mask = bit - 1
    for ket_zero_based in 0:(dimension - 1)
        ket_zero_based & bit == 0 && continue
        bra_zero_based = ket_zero_based ⊻ bit
        fermions_below = count_ones(ket_zero_based & lower_mask)
        sign_value = iseven(fermions_below) ? 1.0 : -1.0
        operator[bra_zero_based + 1, ket_zero_based + 1] = sign_value
    end
    return operator
end

function _thermal_expectation(operator, eigenvectors, boltzmann, partition)
    transformed = adjoint(eigenvectors) * operator * eigenvectors
    return real(sum(boltzmann .* diag(transformed)) / partition)
end

function _merge_peaks(frequencies, weights; tolerance=1.0e-10)
    permutation = sortperm(frequencies)
    merged_frequencies = Float64[]
    merged_weights = Float64[]
    for index in permutation
        frequency = frequencies[index]
        weight = weights[index]
        weight <= 100eps(Float64) && continue
        if !isempty(merged_frequencies) &&
           abs(frequency - last(merged_frequencies)) <= tolerance
            old_weight = last(merged_weights)
            new_weight = old_weight + weight
            merged_frequencies[end] =
                (old_weight * last(merged_frequencies) + weight * frequency) /
                new_weight
            merged_weights[end] = new_weight
        else
            push!(merged_frequencies, frequency)
            push!(merged_weights, weight)
        end
    end
    return (frequencies=merged_frequencies, weights=merged_weights)
end

"""
    solve_two_orbital_atom(U, mu, beta, tau)

Exactly diagonalize one physical cell containing two independent Hubbard
orbitals with the centered interaction
`U*(n_up-1/2)*(n_down-1/2)-mu*(n_up+n_down)`.

The Hilbert space contains only four fermion modes and has dimension 16. The
returned Green function is averaged over the four spin-orbital modes and uses
the same no-leading-minus convention as the DQMC code.
"""
function solve_two_orbital_atom(
    U::Real,
    mu::Real,
    beta::Real,
    tau,
)
    isfinite(U) && U >= 0 || throw(ArgumentError("U must be nonnegative"))
    isfinite(mu) || throw(ArgumentError("mu must be finite"))
    isfinite(beta) && beta > 0 || throw(ArgumentError("beta must be positive"))
    tau_values = Float64.(collect(tau))
    all(time -> 0 <= time <= beta, tau_values) ||
        throw(ArgumentError("tau must lie in [0,beta]"))

    # Mode order: x-up, y-up, x-down, y-down.
    number_of_modes = 4
    annihilation = [
        _annihilation_operator(mode, number_of_modes)
        for mode in 1:number_of_modes
    ]
    creation = adjoint.(annihilation)
    number = [creation[mode] * annihilation[mode] for mode in 1:number_of_modes]
    identity_fock = Matrix{Float64}(I, 16, 16)

    hamiltonian = zeros(Float64, 16, 16)
    for orbital in 1:2
        up = number[orbital]
        down = number[orbital + 2]
        hamiltonian .+= U .* ((up .- 0.5 .* identity_fock) *
                              (down .- 0.5 .* identity_fock))
        hamiltonian .-= mu .* (up + down)
    end
    decomposition = eigen(Hermitian(hamiltonian))
    energies = decomposition.values
    vectors = decomposition.vectors
    energy_shift = minimum(energies)
    boltzmann = exp.(-beta .* (energies .- energy_shift))
    partition_scaled = sum(boltzmann)

    density_operator = sum(number)
    double_occupancy_operator =
        number[1] * number[3] + number[2] * number[4]
    O_alm = (number[1] - number[3] - number[2] + number[4]) / 2
    density = _thermal_expectation(
        density_operator,
        vectors,
        boltzmann,
        partition_scaled,
    )
    double_occupancy = _thermal_expectation(
        double_occupancy_operator,
        vectors,
        boltzmann,
        partition_scaled,
    )
    delta_alm = _thermal_expectation(
        O_alm,
        vectors,
        boltzmann,
        partition_scaled,
    )
    S_alm_q0 = _thermal_expectation(
        O_alm * O_alm,
        vectors,
        boltzmann,
        partition_scaled,
    )

    green = zeros(Float64, length(tau_values))
    peak_frequencies = Float64[]
    peak_weights = Float64[]
    for mode in 1:number_of_modes
        annihilation_eigen = adjoint(vectors) * annihilation[mode] * vectors
        for initial in eachindex(energies), final in eachindex(energies)
            amplitude_squared = abs2(annihilation_eigen[initial, final])
            amplitude_squared <= 100eps(Float64) && continue
            transition = energies[final] - energies[initial]
            for (time_index, time) in pairs(tau_values)
                green[time_index] +=
                    boltzmann[initial] * exp(-time * transition) *
                    amplitude_squared /
                    (partition_scaled * number_of_modes)
            end
            spectral_weight =
                (boltzmann[initial] + boltzmann[final]) *
                amplitude_squared /
                (partition_scaled * number_of_modes)
            push!(peak_frequencies, transition)
            push!(peak_weights, spectral_weight)
        end
    end
    peaks = _merge_peaks(peak_frequencies, peak_weights)
    return (
        U=Float64(U),
        mu=Float64(mu),
        beta=Float64(beta),
        tau=tau_values,
        green=green,
        density=density,
        double_occupancy=double_occupancy,
        delta_alm=delta_alm,
        S_alm_q0=S_alm_q0,
        m_alm_squared=S_alm_q0,
        peak_frequencies=peaks.frequencies,
        peak_weights=peaks.weights,
        spectral_zeroth_moment=sum(peaks.weights),
        hilbert_dimension=16,
    )
end

"""Broaden exact delta peaks by normalized Lorentzians."""
function broaden_spectral_peaks(result, omega; eta::Real=0.1)
    isfinite(eta) && eta > 0 || throw(ArgumentError("eta must be positive"))
    omega_values = Float64.(collect(omega))
    spectrum = zeros(Float64, length(omega_values))
    for (frequency, weight) in zip(
        result.peak_frequencies,
        result.peak_weights,
    )
        @. spectrum += weight * eta /
            (pi * ((omega_values - frequency)^2 + eta^2))
    end
    return (omega=omega_values, spectrum=spectrum, eta=Float64(eta))
end

function save_atom_ed_result(path::AbstractString, result; overwrite=false)
    isfile(path) && !overwrite && throw(ArgumentError(
        "refusing to overwrite existing output: $path",
    ))
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        serialize(io, Dict(
            "schema_version" => 1,
            "kind" => "two_orbital_hubbard_atom_ed",
            "result" => result,
        ))
    end
    return abspath(path)
end

end # module
