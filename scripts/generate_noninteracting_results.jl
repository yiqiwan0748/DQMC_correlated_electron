ENV["GKSwstype"] = "100"

using AltermagnetDQMCTrial
using Plots
using Printf

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT_DIRECTORY = joinpath(REPOSITORY_ROOT, "results", "noninteracting")
const BAND_POINTS_PER_SEGMENT = 160
const CONVERGENCE_NKS = (200, 400, 800)
const PRIMARY_NK = 800
const ENERGY_BINS = 1600
const SMOOTHING_WIDTH = 0.05
const FILLING_POINTS = 2001
const JOURNAL_VALUES = (-0.145, -0.396)

csv_cell(value::Integer) = string(value)
csv_cell(value::AbstractFloat) = @sprintf("%.16g", value)
csv_cell(value) = string(value)

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(csv_cell.(collect(row)), ","))
        end
    end
end

function analyze_grid(Nk)
    samples = sample_band_energies(Nk)
    spectrum = sorted_band_energies(samples)
    dos = density_of_states(
        samples;
        bins=ENERGY_BINS,
        smoothing_width=SMOOTHING_WIDTH,
    )
    peaks = estimate_van_hove_peaks(
        dos,
        samples;
        sorted_energies=spectrum,
    )
    return (samples=samples, spectrum=spectrum, dos=dos, peaks=peaks)
end

mkpath(OUTPUT_DIRECTORY)
band = band_path(points_per_segment=BAND_POINTS_PER_SEGMENT)
analyses = Dict(Nk => analyze_grid(Nk) for Nk in CONVERGENCE_NKS)
primary = analyses[PRIMARY_NK]
inference = infer_doping_convention(
    primary.peaks;
    reported_values=JOURNAL_VALUES,
)

band_rows = (
    (
        band.coordinate[i],
        band.kx[i],
        band.ky[i],
        band.lower_energy[i],
        band.upper_energy[i],
        something(findfirst(==(i), band.tick_indices), 0) == 0 ? "" :
            band.tick_labels[findfirst(==(i), band.tick_indices)],
    ) for i in eachindex(band.coordinate)
)
write_csv(
    joinpath(OUTPUT_DIRECTORY, "band_structure.csv"),
    ("path_coordinate", "kx", "ky", "lower_energy", "upper_energy",
     "high_symmetry_label"),
    band_rows,
)

dos_rows = (
    (PRIMARY_NK, ENERGY_BINS, SMOOTHING_WIDTH, energy, density) for
    (energy, density) in zip(primary.dos.energy, primary.dos.density)
)
write_csv(
    joinpath(OUTPUT_DIRECTORY, "dos.csv"),
    ("Nk", "energy_bins", "smoothing_width", "energy", "dos_per_spin"),
    dos_rows,
)

energy_min, energy_max = extrema(primary.samples.energies)
margin = 0.05 * (energy_max - energy_min)
chemical_potentials = collect(range(
    energy_min - margin,
    energy_max + margin;
    length=FILLING_POINTS,
))
fillings = zero_temperature_filling(
    primary.samples,
    chemical_potentials;
    sorted_energies=primary.spectrum,
)
filling_rows = (
    (PRIMARY_NK, chemical_potential, filling) for
    (chemical_potential, filling) in zip(chemical_potentials, fillings)
)
write_csv(
    joinpath(OUTPUT_DIRECTORY, "filling_curve.csv"),
    ("Nk", "chemical_potential", "filling_including_spin"),
    filling_rows,
)

journal_labels = ("delta_vH_2", "delta_vH")
journal_sorted = sort(collect(JOURNAL_VALUES))
vhs_rows = (
    (
        Nk,
        ENERGY_BINS,
        SMOOTHING_WIDTH,
        peak.order,
        peak.energy,
        peak.filling,
        peak.delta_a,
        peak.delta_b,
        peak.delta_c,
        journal_labels[peak.order],
        journal_sorted[peak.order],
        "delta_c_inferred",
    ) for Nk in CONVERGENCE_NKS for peak in analyses[Nk].peaks
)
write_csv(
    joinpath(OUTPUT_DIRECTORY, "vhs_summary.csv"),
    ("Nk", "energy_bins", "smoothing_width", "peak_order", "mu_vH",
     "filling", "delta_a", "delta_b", "delta_c", "journal_label",
     "journal_value", "inferred_convention"),
    vhs_rows,
)

open(joinpath(OUTPUT_DIRECTORY, "metadata.txt"), "w") do io
    println(io, "milestone=2")
    println(io, "model=two-orbital noninteracting Bloch Hamiltonian")
    println(io, "orbital_basis=(d_xz,d_yz)")
    println(io, "hoppings=t1:-1.0,t2:-1.75,t3:-0.85,t4:-0.65")
    println(io, "momentum_units=radians; lattice_spacing=1")
    println(io, "band_path=Gamma(0,0)-X(pi,0)-M(pi,pi)-Gamma(0,0)")
    println(io, "band_points_per_segment=$BAND_POINTS_PER_SEGMENT")
    println(io, "momentum_grid=k_j=-pi+2pi*j/Nk, j=0,...,Nk-1 (periodic)")
    println(io, "convergence_Nk=$(join(CONVERGENCE_NKS, ','))")
    println(io, "primary_Nk=$PRIMARY_NK")
    println(io, "energy_bins=$ENERGY_BINS")
    println(io, "smoothing_width=$SMOOTHING_WIDTH (Gaussian sigma in energy units)")
    println(io, "smoothing_kernel=truncated at 4 sigma; DOS renormalized after convolution")
    println(io, "DOS=per spin; integral approaches 2 states per physical site")
    println(io, "filling=total electrons per physical site including spin; 0<=n<=4; half filling n=2")
    println(io, "occupied_state_rule=band_energy<=chemical_potential")
    println(io, "vhs_method=two highest separated local DOS maxima after excluding 5 percent of each energy edge; quadratic bin interpolation")
    println(io, "doping_candidates=delta_a:n-2; delta_b:(n-2)/2; delta_c:(n-2)/4")
    println(io, "journal_values=delta_vH:-0.145; delta_vH_2:-0.396")
    println(io, "journal_doping_convention=delta_c (inferred numerically, not independently defined)")
    println(io, "inference_rms_error=$(inference.rms_error)")
end

default(; linewidth=2, framestyle=:box, legend=:best)
band_plot = plot(
    band.coordinate,
    band.lower_energy;
    label="lower band",
    color=:blue,
    xlabel="momentum path",
    ylabel="Energy E/t",
    title="Noninteracting band structure",
    xticks=(band.tick_positions, band.tick_labels),
    size=(900, 550),
    dpi=150,
)
plot!(band_plot, band.coordinate, band.upper_energy; label="upper band", color=:red)
for position in band.tick_positions[2:(end - 1)]
    vline!(band_plot, [position]; color=:gray, alpha=0.4, label=false, linewidth=1)
end
savefig(band_plot, joinpath(OUTPUT_DIRECTORY, "band_structure.png"))

dos_plot = plot(
    primary.dos.energy,
    primary.dos.density;
    label="DOS per spin",
    color=:black,
    xlabel="Energy E/t",
    ylabel="DOS",
    title="DOS (Nk=$PRIMARY_NK)",
)
filling_plot = plot(
    chemical_potentials,
    fillings;
    label="n(mu)",
    color=:darkgreen,
    xlabel="Chemical potential \u03bc/t",
    ylabel="filling n",
    title="Zero-temperature filling (including spin)",
    ylim=(-0.1, 4.1),
)
hline!(filling_plot, [2.0]; linestyle=:dash, color=:gray, label="half filling")
vhs_labels = (
    "second vHs, \u03b4_vH,2 = -0.396",
    "first vHs, \u03b4_vH = -0.145",
)
vhs_colors = ("#0072B2", "#D55E00")
for peak in primary.peaks
    label = vhs_labels[peak.order]
    color = vhs_colors[peak.order]
    vline!(dos_plot, [peak.energy]; linestyle=:dash, color=color, label=label)
    scatter!(
        filling_plot,
        [peak.energy],
        [peak.filling];
        marker=:circle,
        markersize=5,
        color=color,
        label=label,
    )
end
combined_plot = plot(dos_plot, filling_plot; layout=(2, 1), size=(900, 850), dpi=150)
savefig(combined_plot, joinpath(OUTPUT_DIRECTORY, "dos_and_filling.png"))

@printf("band_range_Nk_%d=%.12f,%.12f\n", PRIMARY_NK, energy_min, energy_max)
for Nk in CONVERGENCE_NKS
    for peak in analyses[Nk].peaks
        @printf(
            "Nk=%d peak=%d mu_vH=%.12f filling=%.12f delta_c=%.12f\n",
            Nk,
            peak.order,
            peak.energy,
            peak.filling,
            peak.delta_c,
        )
    end
end
for (coarse_Nk, fine_Nk) in zip(CONVERGENCE_NKS[1:(end - 1)], CONVERGENCE_NKS[2:end])
    for order in 1:2
        coarse = analyses[coarse_Nk].peaks[order]
        fine = analyses[fine_Nk].peaks[order]
        @printf(
            "convergence=%d_to_%d peak=%d abs_delta_energy=%.12g abs_delta_filling=%.12g\n",
            coarse_Nk,
            fine_Nk,
            order,
            abs(fine.energy - coarse.energy),
            abs(fine.filling - coarse.filling),
        )
    end
end
@printf(
    "inferred_convention=%s rms_error=%.12g\n",
    String(inference.convention),
    inference.rms_error,
)
println("output_directory=$OUTPUT_DIRECTORY")
