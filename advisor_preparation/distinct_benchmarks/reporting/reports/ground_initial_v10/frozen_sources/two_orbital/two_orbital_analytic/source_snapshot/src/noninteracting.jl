"""
    ModelParameters{T<:Real}

Immutable hopping parameters for the two-orbital noninteracting model.
"""
struct ModelParameters{T<:Real}
    t1::T
    t2::T
    t3::T
    t4::T
end

function ModelParameters(t1::Real, t2::Real, t3::Real, t4::Real)
    promoted = promote(t1, t2, t3, t4)
    T = typeof(promoted[1])
    return ModelParameters{T}(promoted...)
end

"""
    ModelParameters(; t1=-1.0, t2=-1.75, t3=-0.85, t4=-0.65)

Construct hopping parameters, using the published `Float64` values by default.
"""
ModelParameters(; t1::Real=-1.0, t2::Real=-1.75, t3::Real=-0.85, t4::Real=-0.65) =
    ModelParameters(t1, t2, t3, t4)

const DEFAULT_PARAMETERS = ModelParameters()

"""
    bloch_hamiltonian(kx, ky, parameters=DEFAULT_PARAMETERS)

Return the `2 x 2` real symmetric Bloch Hamiltonian in the orbital basis
`(d_xz, d_yz)`.
"""
function bloch_hamiltonian(
    kx::Real,
    ky::Real,
    parameters::ModelParameters=DEFAULT_PARAMETERS,
)
    (; t1, t2, t3, t4) = parameters
    cos_kx = cos(kx)
    cos_ky = cos(ky)
    sin_kx = sin(kx)
    sin_ky = sin(ky)

    epsilon_x = -2 * t1 * cos_kx - 2 * t2 * cos_ky - 4 * t3 * cos_kx * cos_ky
    epsilon_y = -2 * t2 * cos_kx - 2 * t1 * cos_ky - 4 * t3 * cos_kx * cos_ky
    epsilon_xy = -4 * t4 * sin_kx * sin_ky

    return [epsilon_x epsilon_xy; epsilon_xy epsilon_y]
end

"""
    band_energies(kx, ky, parameters=DEFAULT_PARAMETERS)

Return the two eigenvalues of the Bloch Hamiltonian in ascending order.
"""
function band_energies(
    kx::Real,
    ky::Real,
    parameters::ModelParameters=DEFAULT_PARAMETERS,
)
    hamiltonian = bloch_hamiltonian(kx, ky, parameters)
    epsilon_x = hamiltonian[1, 1]
    epsilon_y = hamiltonian[2, 2]
    epsilon_xy = hamiltonian[1, 2]

    center = epsilon_x / 2 + epsilon_y / 2
    radius = hypot(epsilon_x / 2 - epsilon_y / 2, epsilon_xy)
    return [center - radius, center + radius]
end

const DEFAULT_BAND_PATH_VERTICES = (
    (0.0, 0.0),
    (pi, 0.0),
    (pi, pi),
    (0.0, 0.0),
)
const DEFAULT_BAND_PATH_LABELS = ("Γ", "X", "M", "Γ")

"""
    band_path(; vertices=DEFAULT_BAND_PATH_VERTICES,
                labels=DEFAULT_BAND_PATH_LABELS,
                points_per_segment=100,
                parameters=DEFAULT_PARAMETERS)

Sample both bands on a piecewise-linear path (default `Γ-X-M-Γ`). Shared
segment endpoints occur only once. The result includes coordinates, momenta,
energies, and high-symmetry tick positions, labels, and indices.
"""
function band_path(;
    vertices=DEFAULT_BAND_PATH_VERTICES,
    labels=DEFAULT_BAND_PATH_LABELS,
    points_per_segment=100,
    parameters::ModelParameters=DEFAULT_PARAMETERS,
)
    path_vertices = collect(vertices)
    path_labels = String.(collect(labels))
    length(path_vertices) >= 2 ||
        throw(ArgumentError("a band path needs at least two vertices"))
    length(path_labels) == length(path_vertices) ||
        throw(ArgumentError("labels and vertices must have the same length"))
    all(vertex -> length(vertex) == 2, path_vertices) ||
        throw(ArgumentError("each path vertex must contain kx and ky"))

    segment_count = length(path_vertices) - 1
    samples_per_segment = if points_per_segment isa Integer
        fill(Int(points_per_segment), segment_count)
    else
        Int.(collect(points_per_segment))
    end
    length(samples_per_segment) == segment_count || throw(
        ArgumentError("points_per_segment must have one entry per path segment"),
    )
    all(>(0), samples_per_segment) ||
        throw(ArgumentError("points_per_segment entries must be positive"))

    coordinate_types = [
        typeof(float(component)) for vertex in path_vertices for component in vertex
    ]
    T = promote_type(typeof(float(parameters.t1)), coordinate_types...)
    converted_vertices = [
        (convert(T, vertex[1]), convert(T, vertex[2])) for vertex in path_vertices
    ]
    segment_lengths = [
        hypot(
            converted_vertices[i + 1][1] - converted_vertices[i][1],
            converted_vertices[i + 1][2] - converted_vertices[i][2],
        ) for i in 1:segment_count
    ]
    tick_positions = Vector{T}(undef, length(converted_vertices))
    tick_positions[1] = zero(T)
    for i in 1:segment_count
        tick_positions[i + 1] = tick_positions[i] + segment_lengths[i]
    end
    tick_indices = vcat(1, 1 .+ cumsum(samples_per_segment))

    total_points = 1 + sum(samples_per_segment)
    coordinate = Vector{T}(undef, total_points)
    kx = Vector{T}(undef, total_points)
    ky = Vector{T}(undef, total_points)
    first_energies = band_energies(converted_vertices[1]..., parameters)
    E = eltype(first_energies)
    lower_energy = Vector{E}(undef, total_points)
    upper_energy = Vector{E}(undef, total_points)

    coordinate[1] = zero(T)
    kx[1], ky[1] = converted_vertices[1]
    lower_energy[1], upper_energy[1] = first_energies
    output_index = 1
    for segment in 1:segment_count
        start_vertex = converted_vertices[segment]
        end_vertex = converted_vertices[segment + 1]
        segment_samples = samples_per_segment[segment]
        for j in 1:segment_samples
            output_index += 1
            if j == segment_samples
                kx[output_index], ky[output_index] = end_vertex
                coordinate[output_index] = tick_positions[segment + 1]
            else
                fraction = convert(T, j) / convert(T, segment_samples)
                kx[output_index] =
                    start_vertex[1] + fraction * (end_vertex[1] - start_vertex[1])
                ky[output_index] =
                    start_vertex[2] + fraction * (end_vertex[2] - start_vertex[2])
                coordinate[output_index] =
                    tick_positions[segment] + fraction * segment_lengths[segment]
            end
            energies = band_energies(kx[output_index], ky[output_index], parameters)
            lower_energy[output_index], upper_energy[output_index] = energies
        end
    end

    return (
        coordinate=coordinate,
        kx=kx,
        ky=ky,
        lower_energy=lower_energy,
        upper_energy=upper_energy,
        tick_positions=tick_positions,
        tick_labels=path_labels,
        tick_indices=tick_indices,
    )
end

"""
Band energies on a periodic square grid. `k_values` are `-π + 2πj/Nk` for
`j = 0, ..., Nk-1`. `energies` has lower/upper rows and `Nk^2` columns, with
`kx` varying fastest.
"""
struct BandGridSamples{K<:AbstractFloat,E<:Real}
    Nk::Int
    k_values::Vector{K}
    energies::Matrix{E}
end

"""
    sample_band_energies(Nk; parameters=DEFAULT_PARAMETERS)

Evaluate both bands once on `[-π, π) × [-π, π)`. Reuse the result for DOS,
filling, and van-Hove calculations.
"""
function sample_band_energies(
    Nk::Integer;
    parameters::ModelParameters=DEFAULT_PARAMETERS,
)
    Nk > 0 || throw(ArgumentError("Nk must be positive"))
    grid_size = Int(Nk)
    K = typeof(float(parameters.t1))
    pi_value = convert(K, pi)
    momentum_step = 2 * pi_value / convert(K, grid_size)
    k_values = [
        -pi_value + convert(K, j) * momentum_step for j in 0:(grid_size - 1)
    ]
    first_energies = band_energies(k_values[1], k_values[1], parameters)
    energies = Matrix{eltype(first_energies)}(undef, 2, grid_size^2)
    for iy in 1:grid_size, ix in 1:grid_size
        column = ix + (iy - 1) * grid_size
        energies[:, column] = band_energies(k_values[ix], k_values[iy], parameters)
    end
    return BandGridSamples(grid_size, k_values, energies)
end

sample_band_energies(Nk::Integer, parameters::ModelParameters) =
    sample_band_energies(Nk; parameters=parameters)

"""Return a sorted copy of all `2Nk^2` one-spin band energies."""
sorted_band_energies(samples::BandGridSamples) =
    sort!(vec(copy(samples.energies)))

function _gaussian_smooth(values::AbstractVector, bin_width::Real, sigma::Real)
    sigma > 0 || return copy(values)
    radius = min(length(values) - 1, ceil(Int, 4 * sigma / bin_width))
    radius > 0 || return copy(values)
    offsets = collect(-radius:radius)
    kernel = exp.(-0.5 .* ((offsets .* bin_width) ./ sigma) .^ 2)
    kernel ./= sum(kernel)
    smoothed = zeros(promote_type(eltype(values), eltype(kernel)), length(values))
    for i in eachindex(values)
        for (kernel_index, offset) in enumerate(offsets)
            source_index = i - offset
            if firstindex(values) <= source_index <= lastindex(values)
                smoothed[i] += kernel[kernel_index] * values[source_index]
            end
        end
    end
    return smoothed
end

"""
    density_of_states(samples; bins=1000, smoothing_width=0)

Compute a histogram DOS per spin. `smoothing_width` is the standard deviation
of an optional Gaussian convolution in energy units, truncated at four standard
deviations. The midpoint-rule integral is normalized to two states per site.
"""
function density_of_states(
    samples::BandGridSamples;
    bins::Integer=1000,
    smoothing_width::Real=0,
)
    bins >= 2 || throw(ArgumentError("bins must be at least two"))
    isfinite(smoothing_width) && smoothing_width >= 0 ||
        throw(ArgumentError("smoothing_width must be finite and nonnegative"))
    bin_count = Int(bins)
    energy_min, energy_max = extrema(samples.energies)
    energy_max > energy_min ||
        throw(ArgumentError("the sampled spectrum must have nonzero width"))
    bin_width = (energy_max - energy_min) / bin_count
    bin_edges = collect(range(energy_min, energy_max; length=bin_count + 1))
    energy = [energy_min + (i - 0.5) * bin_width for i in 1:bin_count]
    counts = zeros(Int, bin_count)
    for sampled_energy in samples.energies
        bin_index = sampled_energy == energy_max ? bin_count :
                    floor(Int, (sampled_energy - energy_min) / bin_width) + 1
        counts[clamp(bin_index, 1, bin_count)] += 1
    end
    density = counts ./ (samples.Nk^2 * bin_width)
    density = _gaussian_smooth(density, bin_width, smoothing_width)
    # Correct the small edge loss from a truncated smoothing convolution.
    density .*= 2 / (sum(density) * bin_width)
    return (
        energy=energy,
        density=density,
        counts=counts,
        bin_edges=bin_edges,
        bin_width=bin_width,
        smoothing_width=smoothing_width,
        Nk=samples.Nk,
    )
end

function _checked_sorted_energies(samples::BandGridSamples, spectrum)
    if isnothing(spectrum)
        return sorted_band_energies(samples)
    end
    length(spectrum) == length(samples.energies) ||
        throw(ArgumentError("sorted_energies has the wrong length"))
    issorted(spectrum) || throw(ArgumentError("sorted_energies must be sorted"))
    return spectrum
end

"""
    zero_temperature_filling(samples, chemical_potential;
                             sorted_energies=nothing)

Return electron filling per physical site, including spin degeneracy. Occupied
states satisfy `energy <= chemical_potential`; the range is zero to four and
half filling is two. Supply one precomputed sorted spectrum for repeated calls.
"""
function zero_temperature_filling(
    samples::BandGridSamples,
    chemical_potential::Real;
    sorted_energies=nothing,
)
    spectrum = _checked_sorted_energies(samples, sorted_energies)
    return 2 * searchsortedlast(spectrum, chemical_potential) / samples.Nk^2
end

function zero_temperature_filling(
    samples::BandGridSamples,
    chemical_potentials::AbstractVector{<:Real};
    sorted_energies=nothing,
)
    spectrum = _checked_sorted_energies(samples, sorted_energies)
    normalization = 2 / samples.Nk^2
    return [
        normalization * searchsortedlast(spectrum, chemical_potential) for
        chemical_potential in chemical_potentials
    ]
end

function _quadratic_peak(energy, density, index)
    left = density[index - 1]
    center = density[index]
    right = density[index + 1]
    curvature = left - 2 * center + right
    offset = iszero(curvature) ? zero(curvature) :
             0.5 * (left - right) / curvature
    offset = clamp(offset, -0.5, 0.5)
    spacing = (energy[index + 1] - energy[index - 1]) / 2
    peak_energy = energy[index] + offset * spacing
    peak_density = center +
                   0.5 * (right - left) * offset +
                   0.5 * curvature * offset^2
    return peak_energy, peak_density
end

"""
    estimate_van_hove_peaks(dos, samples; number=2, edge_fraction=0.05,
                            min_separation=nothing, sorted_energies=nothing)

Estimate principal internal DOS maxima. Candidates are ranked by height after
excluding band edges and enforcing energy separation. Three-point quadratic
interpolation refines each maximum. Results include `mu_vH`, filling, and the
three requested candidate doping conventions.
"""
function estimate_van_hove_peaks(
    dos,
    samples::BandGridSamples;
    number::Integer=2,
    edge_fraction::Real=0.05,
    min_separation=nothing,
    sorted_energies=nothing,
)
    number > 0 || throw(ArgumentError("number must be positive"))
    0 <= edge_fraction < 0.5 ||
        throw(ArgumentError("edge_fraction must lie in [0, 0.5)"))
    energy = dos.energy
    density = dos.density
    length(energy) == length(density) >= 3 ||
        throw(ArgumentError("DOS energy and density arrays must match"))
    energy_range = last(energy) - first(energy)
    separation = isnothing(min_separation) ? 0.05 * energy_range : min_separation
    separation >= 0 || throw(ArgumentError("min_separation must be nonnegative"))
    lower_limit = first(energy) + edge_fraction * energy_range
    upper_limit = last(energy) - edge_fraction * energy_range

    candidates = [
        i for i in 2:(length(density) - 1) if
        lower_limit <= energy[i] <= upper_limit &&
        density[i] >= density[i - 1] && density[i] > density[i + 1]
    ]
    sort!(candidates; by=i -> (-density[i], energy[i]))
    selected = Int[]
    for candidate in candidates
        separated = all(
            abs(energy[candidate] - energy[chosen]) >= separation for
            chosen in selected
        )
        if separated
            push!(selected, candidate)
            length(selected) == number && break
        end
    end
    length(selected) == number || throw(
        ArgumentError("fewer than $number separated internal DOS peaks were found"),
    )
    sort!(selected; by=i -> energy[i])

    spectrum = _checked_sorted_energies(samples, sorted_energies)
    normalization = 2 / samples.Nk^2
    peaks = NamedTuple[]
    for (order, index) in enumerate(selected)
        peak_energy, peak_density = _quadratic_peak(energy, density, index)
        filling = normalization * searchsortedlast(spectrum, peak_energy)
        push!(peaks, (
            order=order,
            bin_index=index,
            energy=peak_energy,
            density=peak_density,
            filling=filling,
            delta_a=filling - 2,
            delta_b=(filling - 2) / 2,
            delta_c=(filling - 2) / 4,
        ))
    end
    return peaks
end

"""
    infer_doping_convention(peaks; reported_values=(-0.145, -0.396))

Compare sorted candidate dopings with the journal values. The smallest RMS
mismatch identifies a numerical inference, not an independent definition of
the journal's notation.
"""
function infer_doping_convention(
    peaks;
    reported_values=(-0.145, -0.396),
)
    length(peaks) == length(reported_values) || throw(
        ArgumentError("the peak and reported-value counts must match"),
    )
    reported = sort(float.(collect(reported_values)))
    function mismatch(field)
        estimated = sort([float(getproperty(peak, field)) for peak in peaks])
        return sqrt(sum(abs2, estimated .- reported) / length(reported))
    end
    errors = (
        delta_a=mismatch(:delta_a),
        delta_b=mismatch(:delta_b),
        delta_c=mismatch(:delta_c),
    )
    conventions = (:delta_a, :delta_b, :delta_c)
    best = conventions[argmin([
        getproperty(errors, convention) for convention in conventions
    ])]
    return (convention=best, rms_error=getproperty(errors, best), errors=errors)
end
