ENV["GKSwstype"] = "100"

using AltermagnetDQMCTrial
using Plots
using Printf

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT_DIRECTORY = joinpath(
    REPOSITORY_ROOT,
    "results",
    "moreo_L100_noninteracting",
)

const BAND_FIGURE_PATH = joinpath(
    OUTPUT_DIRECTORY,
    "bands_Gamma_X_M_Gamma_L100.png",
)
const FERMI_SURFACE_FIGURE_PATH = joinpath(
    OUTPUT_DIRECTORY,
    "fermi_surface_L100.png",
)
const COMBINED_FIGURE_PATH = joinpath(
    OUTPUT_DIRECTORY,
    "bands_and_fermi_surface_L100.png",
)
const COMPARISON_PATH = joinpath(
    OUTPUT_DIRECTORY,
    "finite_size_comparison_L4_L100.csv",
)
const BAND_PATH_DATA_PATH = joinpath(OUTPUT_DIRECTORY, "band_path.csv")
const METADATA_PATH = joinpath(OUTPUT_DIRECTORY, "metadata.txt")

const L_SMALL = 4
const L_LARGE = 100
const BETA = 4.0
const MU = -1.897
const INTERACTION_U = 0.0
const PARAMETERS = ModelParameters()
const SPIN_DEGENERACY = 2

const BAND_POINTS_PER_SEGMENT = 400
const VISUALIZATION_MESH_POINTS = 601
const VALIDATION_TOLERANCE = 1.0e-12

const LOWER_BAND_COLOR = "#0072B2"
const UPPER_BAND_COLOR = "#D55E00"
const INK_COLOR = "#1F2937"
const GUIDE_COLOR = "#94A3B8"

const L4_REFERENCE = (
    density=1.820507020579327,
    kinetic_energy_per_cell=-5.402514621495678,
    grand_energy_per_cell=-1.949012803456695,
)

const HIGH_SYMMETRY_REFERENCES = (
    (label="Gamma", kx=0.0, ky=0.0, lower=8.9, upper=8.9),
    (label="X", kx=pi, ky=0.0, lower=-4.9, upper=-1.9),
    (label="M", kx=pi, ky=pi, lower=-2.1, upper=-2.1),
)

function stable_fermi(scaled_energy::T) where {T<:AbstractFloat}
    if scaled_energy >= zero(T)
        tail = exp(-scaled_energy)
        return tail / (one(T) + tail)
    else
        tail = exp(scaled_energy)
        return one(T) / (one(T) + tail)
    end
end

function audit_discrete_grid(L)
    occupation_sum = 0.0
    kinetic_sum = 0.0
    grand_sum = 0.0
    minimum_occupation = Inf
    maximum_occupation = -Inf
    all_energies_finite = true
    all_occupations_finite = true
    all_occupations_in_range = true

    for ky in allowed_momenta(L), kx in allowed_momenta(L)
        for energy in band_energies(kx, ky, PARAMETERS)
            occupation = stable_fermi(BETA * (energy - MU))
            all_energies_finite &= isfinite(energy)
            all_occupations_finite &= isfinite(occupation)
            all_occupations_in_range &= 0.0 <= occupation <= 1.0
            minimum_occupation = min(minimum_occupation, occupation)
            maximum_occupation = max(maximum_occupation, occupation)
            occupation_sum += occupation
            kinetic_sum += energy * occupation
            grand_sum += (energy - MU) * occupation
        end
    end

    normalization = SPIN_DEGENERACY / L^2
    return (
        density=normalization * occupation_sum,
        kinetic_energy_per_cell=normalization * kinetic_sum,
        grand_energy_per_cell=normalization * grand_sum,
        minimum_occupation=minimum_occupation,
        maximum_occupation=maximum_occupation,
        all_energies_finite=all_energies_finite,
        all_occupations_finite=all_occupations_finite,
        all_occupations_in_range=all_occupations_in_range,
        momentum_count=L^2,
        one_spin_level_count=2 * L^2,
    )
end

function check_thermodynamics(L)
    thermodynamics = finite_lattice_free_thermodynamics(
        L,
        BETA,
        MU;
        parameters=PARAMETERS,
    )
    audit = audit_discrete_grid(L)
    observable_values = (
        thermodynamics.density,
        thermodynamics.kinetic_energy_per_cell,
        thermodynamics.grand_energy_per_cell,
    )
    all_observables_finite = all(isfinite, observable_values)
    density_in_range = 0.0 <= thermodynamics.density <= 4.0
    identity_error = abs(
        thermodynamics.grand_energy_per_cell -
        (thermodynamics.kinetic_energy_per_cell - MU * thermodynamics.density),
    )
    formula_errors = (
        density=abs(thermodynamics.density - audit.density),
        kinetic_energy_per_cell=abs(
            thermodynamics.kinetic_energy_per_cell -
            audit.kinetic_energy_per_cell,
        ),
        grand_energy_per_cell=abs(
            thermodynamics.grand_energy_per_cell -
            audit.grand_energy_per_cell,
        ),
    )

    audit.all_energies_finite || error("L=$L band energy is not finite")
    audit.all_occupations_finite || error("L=$L occupation is not finite")
    audit.all_occupations_in_range || error("L=$L occupation lies outside [0,1]")
    all_observables_finite || error("L=$L thermodynamic observable is not finite")
    density_in_range || error("L=$L density lies outside [0,4]")
    identity_error <= VALIDATION_TOLERANCE || error(
        "L=$L thermodynamic identity failed with error $identity_error",
    )
    maximum(Tuple(formula_errors)) <= VALIDATION_TOLERANCE || error(
        "L=$L validated API and direct formula audit disagree",
    )

    return (
        thermodynamics=thermodynamics,
        audit=audit,
        all_observables_finite=all_observables_finite,
        density_in_range=density_in_range,
        identity_error=identity_error,
        formula_errors=formula_errors,
    )
end

function check_l4_reference(thermodynamics)
    errors = (
        density=abs(thermodynamics.density - L4_REFERENCE.density),
        kinetic_energy_per_cell=abs(
            thermodynamics.kinetic_energy_per_cell -
            L4_REFERENCE.kinetic_energy_per_cell,
        ),
        grand_energy_per_cell=abs(
            thermodynamics.grand_energy_per_cell -
            L4_REFERENCE.grand_energy_per_cell,
        ),
    )
    maximum(Tuple(errors)) <= VALIDATION_TOLERANCE || error(
        "new L=4 thermodynamics do not reproduce the validated reference",
    )
    return errors
end

function check_high_symmetry_points()
    rows = map(HIGH_SYMMETRY_REFERENCES) do reference
        energies = band_energies(reference.kx, reference.ky, PARAMETERS)
        error = max(
            abs(energies[1] - reference.lower),
            abs(energies[2] - reference.upper),
        )
        return merge(
            reference,
            (
                actual_lower=energies[1],
                actual_upper=energies[2],
                maximum_abs_error=error,
            ),
        )
    end
    maximum_error = maximum(row.maximum_abs_error for row in rows)
    maximum_error <= VALIDATION_TOLERANCE || error(
        "Milestone 1 high-symmetry values disagree by $maximum_error",
    )
    return rows, maximum_error
end

function check_band_path_high_symmetry(path)
    path_references = (
        HIGH_SYMMETRY_REFERENCES...,
        HIGH_SYMMETRY_REFERENCES[1],
    )
    errors = [
        max(
            abs(path.lower_energy[index] - reference.lower),
            abs(path.upper_energy[index] - reference.upper),
        ) for (index, reference) in zip(path.tick_indices, path_references)
    ]
    maximum_error = maximum(errors)
    maximum_error <= VALIDATION_TOLERANCE || error(
        "band-path high-symmetry values disagree by $maximum_error",
    )
    return maximum_error
end

function comparison_rows(L4, L100)
    values = (
        (
            observable="electron_density_n",
            L4=L4.density,
            L100=L100.density,
        ),
        (
            observable="kinetic_energy_H0_per_L2",
            L4=L4.kinetic_energy_per_cell,
            L100=L100.kinetic_energy_per_cell,
        ),
        (
            observable="grand_energy_H_minus_muN_per_L2",
            L4=L4.grand_energy_per_cell,
            L100=L100.grand_energy_per_cell,
        ),
    )
    return map(values) do row
        difference = row.L100 - row.L4
        merge(
            row,
            (
                signed_difference_L100_minus_L4=difference,
                absolute_difference=abs(difference),
            ),
        )
    end
end

function write_comparison(path, rows)
    open(path, "w") do io
        println(
            io,
            "observable,L4,L100,signed_difference_L100_minus_L4," *
            "absolute_difference",
        )
        for row in rows
            @printf(
                io,
                "%s,%.17g,%.17g,%.17g,%.17g\n",
                row.observable,
                row.L4,
                row.L100,
                row.signed_difference_L100_minus_L4,
                row.absolute_difference,
            )
        end
    end
end

function write_band_path(path, band)
    symmetry_labels = Dict(
        index => label for (index, label) in
        zip(band.tick_indices, band.tick_labels)
    )
    open(path, "w") do io
        println(
            io,
            "path_index,path_distance,kx,ky,epsilon_minus,epsilon_plus," *
            "symmetry_label",
        )
        for index in eachindex(band.coordinate)
            @printf(
                io,
                "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%s\n",
                index,
                band.coordinate[index],
                band.kx[index],
                band.ky[index],
                band.lower_energy[index],
                band.upper_energy[index],
                get(symmetry_labels, index, ""),
            )
        end
    end
end

function visualization_mesh()
    momenta = collect(range(-pi, pi; length=VISUALIZATION_MESH_POINTS))
    lower = Matrix{Float64}(
        undef,
        VISUALIZATION_MESH_POINTS,
        VISUALIZATION_MESH_POINTS,
    )
    upper = similar(lower)
    for (iy, ky) in enumerate(momenta), (ix, kx) in enumerate(momenta)
        energies = band_energies(kx, ky, PARAMETERS)
        lower[iy, ix] = energies[1]
        upper[iy, ix] = energies[2]
    end
    all(isfinite, lower) || error("lower-band visualization mesh is not finite")
    all(isfinite, upper) || error("upper-band visualization mesh is not finite")
    return momenta, lower, upper
end

function make_band_plot(path)
    energy_min = min(minimum(path.lower_energy), minimum(path.upper_energy))
    energy_max = max(maximum(path.lower_energy), maximum(path.upper_energy))
    padding = 0.04 * (energy_max - energy_min)
    plot_object = plot(
        path.coordinate,
        path.lower_energy;
        label="lower band  ε₋",
        color=LOWER_BAND_COLOR,
        linewidth=3.0,
        xlabel="Wave-vector path",
        ylabel="Band energy ε/t",
        title="Two energy bands along Γ-X-M-Γ",
        xticks=(path.tick_positions, path.tick_labels),
        xlims=(first(path.coordinate), last(path.coordinate)),
        ylims=(energy_min - padding, energy_max + padding),
        legend=:top,
        left_margin=11Plots.mm,
    )
    plot!(
        plot_object,
        path.coordinate,
        path.upper_energy;
        label="upper band  ε₊",
        color=UPPER_BAND_COLOR,
        linewidth=3.0,
    )
    for position in path.tick_positions[2:(end - 1)]
        vline!(
            plot_object,
            [position];
            color=GUIDE_COLOR,
            alpha=0.55,
            linewidth=1.2,
            label=false,
        )
    end
    hline!(
        plot_object,
        [MU];
        label=@sprintf("μ/t = %.3f", MU),
        color=INK_COLOR,
        linestyle=:dash,
        linewidth=2.4,
    )
    return plot_object
end

function make_fermi_surface_plot(momenta, lower, upper)
    lower_gradient = cgrad([LOWER_BAND_COLOR, LOWER_BAND_COLOR])
    upper_gradient = cgrad([UPPER_BAND_COLOR, UPPER_BAND_COLOR])
    plot_object = contour(
        momenta,
        momenta,
        lower;
        levels=[MU],
        color=lower_gradient,
        colorbar=false,
        linewidth=3.2,
        label=false,
        xlabel="kₓ",
        ylabel="kᵧ",
        title="Fermi contours at μ/t = -1.897",
        xlims=(-pi, pi),
        ylims=(-pi, pi),
        xticks=([-pi, 0.0, pi], ["−π", "0", "π"]),
        yticks=([-pi, 0.0, pi], ["−π", "0", "π"]),
        aspect_ratio=:equal,
        legend=(0.5, 0.5),
    )
    contour!(
        plot_object,
        momenta,
        momenta,
        upper;
        levels=[MU],
        color=upper_gradient,
        colorbar=false,
        linewidth=3.2,
        label=false,
    )
    plot!(
        plot_object,
        [NaN],
        [NaN];
        label="lower band: ε₋(k) = μ",
        color=LOWER_BAND_COLOR,
        linewidth=3.2,
    )
    plot!(
        plot_object,
        [NaN],
        [NaN];
        label="upper band: ε₊(k) = μ",
        color=UPPER_BAND_COLOR,
        linewidth=3.2,
    )
    return plot_object
end

function configure_plot_defaults()
    default(
        background_color=:white,
        foreground_color_subplot=INK_COLOR,
        framestyle=:box,
        grid=false,
        guidefontsize=14,
        tickfontsize=12,
        titlefontsize=15,
        legendfontsize=11,
        linewidth=2.4,
        fontfamily="sans-serif",
        margin=7Plots.mm,
    )
end

function make_footer_plot()
    plot_object = plot(
        ;
        xlims=(0.0, 1.0),
        ylims=(0.0, 1.0),
        framestyle=:none,
        showaxis=false,
        ticks=nothing,
        grid=false,
        legend=false,
        margin=0Plots.mm,
    )
    annotate!(
        plot_object,
        0.5,
        0.5,
        text(
            "Contours use a dense mesh for visualization; finite-L thermodynamics use the discrete 100×100 PBC momentum grid at βt=4.0.",
            10,
            INK_COLOR,
            :center,
        ),
    )
    return plot_object
end

function write_metadata(
    path,
    L4,
    L100,
    L4_reference_errors,
    high_symmetry_rows,
    high_symmetry_error,
    band_path_high_symmetry_error,
    generation_mode,
)
    open(path, "w") do io
        println(io, "artifact=Professor Moreo exact noninteracting L=100 reference")
        println(io, "generator=scripts/generate_moreo_L100_noninteracting.jl")
        println(io, "generation_mode=$generation_mode")
        println(io, "calculation=exact finite-grid free fermions; no sampling")
        println(io, "U=$INTERACTION_U")
        println(io, "L_values=$L_SMALL,$L_LARGE")
        println(io, "beta_t=$BETA")
        println(io, "mu_over_t=$MU")
        println(io, "t1=$(PARAMETERS.t1)")
        println(io, "t2=$(PARAMETERS.t2)")
        println(io, "t3=$(PARAMETERS.t3)")
        println(io, "t4=$(PARAMETERS.t4)")
        println(io, "orbital_basis=(d_xz,d_yz)")
        println(io, "lattice=square; lattice_spacing=1")
        println(io, "momentum_units=radians")
        println(io, "energy_units=t")
        println(io, "boundary_conditions=periodic in x and y")
        println(io, "allowed_momenta=k_alpha=2*pi*n_alpha/L, n_alpha=0,...,L-1")
        println(io, "spin_degeneracy=$SPIN_DEGENERACY")
        println(io, "density_normalization=electrons per L^2 unit cells; 0<=n<=4")
        println(io, "band_energy_convention=bare eigenvalues of H0(k); mu is not absorbed")
        println(io, "fermi_function=f_beta(x)=1/(1+exp(beta*x)), x=epsilon-mu")
        println(io, "density_formula=n=(2/L^2)*sum_(k,nu) f_beta(epsilon_nu(k)-mu)")
        println(io, "kinetic_formula=<H0>/L^2=(2/L^2)*sum_(k,nu) epsilon_nu(k)*f_beta")
        println(io, "grand_energy_formula=<H-mu*N>/L^2=(2/L^2)*sum_(k,nu) (epsilon_nu(k)-mu)*f_beta")
        println(io, "grand_energy_interpretation=expectation of H-mu*N; not the grand potential")
        println(io, "finite_l_api=finite_lattice_free_thermodynamics")
        println(io, "band_api=band_energies and band_path")
        println(io, "band_path=Gamma(0,0)-X(pi,0)-M(pi,pi)-Gamma(0,0)")
        println(io, "band_path_coordinate=standard cumulative Euclidean path distance")
        println(io, "band_points_per_segment=$BAND_POINTS_PER_SEGMENT")
        println(io, "band_path_lattice_size_dependence=none")
        println(io, "visualization_mesh=$(VISUALIZATION_MESH_POINTS)x$(VISUALIZATION_MESH_POINTS) inclusive grid on [-pi,pi]^2")
        println(io, "visualization_mesh_role=smooth contour rendering only; not used for thermodynamics")
        println(io, "continuous_fermi_contour_lattice_size_dependence=none")
        println(io, "L100_thermodynamic_grid=100x100 discrete periodic allowed momenta")
        println(io, "dqmc_state_constructed=false")
        println(io, "metropolis_sampling_performed=false")
        println(io, "validation_tolerance=$VALIDATION_TOLERANCE (absolute)")
        println(io, "L4_observables_finite=$(L4.all_observables_finite)")
        println(io, "L100_observables_finite=$(L100.all_observables_finite)")
        println(io, "L4_density_in_range=$(L4.density_in_range)")
        println(io, "L100_density_in_range=$(L100.density_in_range)")
        println(io, "L4_all_band_energies_finite=$(L4.audit.all_energies_finite)")
        println(io, "L100_all_band_energies_finite=$(L100.audit.all_energies_finite)")
        println(io, "L4_all_occupations_finite=$(L4.audit.all_occupations_finite)")
        println(io, "L100_all_occupations_finite=$(L100.audit.all_occupations_finite)")
        println(io, "L4_all_occupations_in_range=$(L4.audit.all_occupations_in_range)")
        println(io, "L100_all_occupations_in_range=$(L100.audit.all_occupations_in_range)")
        @printf(io, "L4_minimum_occupation=%.17g\n", L4.audit.minimum_occupation)
        @printf(io, "L4_maximum_occupation=%.17g\n", L4.audit.maximum_occupation)
        @printf(io, "L100_minimum_occupation=%.17g\n", L100.audit.minimum_occupation)
        @printf(io, "L100_maximum_occupation=%.17g\n", L100.audit.maximum_occupation)
        @printf(io, "L4_identity_error=%.17g\n", L4.identity_error)
        @printf(io, "L100_identity_error=%.17g\n", L100.identity_error)
        @printf(io, "L4_formula_density_error=%.17g\n", L4.formula_errors.density)
        @printf(io, "L4_formula_kinetic_error=%.17g\n", L4.formula_errors.kinetic_energy_per_cell)
        @printf(io, "L4_formula_grand_energy_error=%.17g\n", L4.formula_errors.grand_energy_per_cell)
        @printf(io, "L100_formula_density_error=%.17g\n", L100.formula_errors.density)
        @printf(io, "L100_formula_kinetic_error=%.17g\n", L100.formula_errors.kinetic_energy_per_cell)
        @printf(io, "L100_formula_grand_energy_error=%.17g\n", L100.formula_errors.grand_energy_per_cell)
        println(io, "L4_reference_source=test/free_fermion_tests.jl")
        @printf(io, "L4_density_reference=%.17g\n", L4_REFERENCE.density)
        @printf(io, "L4_kinetic_reference=%.17g\n", L4_REFERENCE.kinetic_energy_per_cell)
        @printf(io, "L4_grand_energy_reference=%.17g\n", L4_REFERENCE.grand_energy_per_cell)
        @printf(io, "L4_density_reference_error=%.17g\n", L4_reference_errors.density)
        @printf(io, "L4_kinetic_reference_error=%.17g\n", L4_reference_errors.kinetic_energy_per_cell)
        @printf(io, "L4_grand_energy_reference_error=%.17g\n", L4_reference_errors.grand_energy_per_cell)
        @printf(io, "high_symmetry_max_abs_error=%.17g\n", high_symmetry_error)
        @printf(io, "band_path_high_symmetry_max_abs_error=%.17g\n", band_path_high_symmetry_error)
        for row in high_symmetry_rows
            @printf(
                io,
                "high_symmetry_%s=epsilon_minus:%.17g,epsilon_plus:%.17g,error:%.17g\n",
                row.label,
                row.actual_lower,
                row.actual_upper,
                row.maximum_abs_error,
            )
        end
        println(io, "output_bands=bands_Gamma_X_M_Gamma_L100.png")
        println(io, "output_fermi_surface=fermi_surface_L100.png")
        println(io, "output_combined=bands_and_fermi_surface_L100.png")
        println(io, "output_comparison=finite_size_comparison_L4_L100.csv")
        println(io, "output_band_path=band_path.csv")
    end
end

function prepare_output_directory(arguments)
    generation_mode = if isempty(arguments)
        "initial"
    elseif arguments == ["--regenerate"]
        "deterministic_regeneration"
    else
        error("usage: julia --project=. scripts/generate_moreo_L100_noninteracting.jl [--regenerate]")
    end

    if generation_mode == "initial"
        !ispath(OUTPUT_DIRECTORY) || error(
            "output directory already exists: $OUTPUT_DIRECTORY",
        )
    else
        isdir(OUTPUT_DIRECTORY) || error(
            "--regenerate requires the existing output directory: $OUTPUT_DIRECTORY",
        )
    end
    return generation_mode
end

function print_report(rows, L4, L100, L4_reference_errors, high_symmetry_error, band_path_error)
    println("comparison_table")
    println("observable,L4,L100,signed_difference_L100_minus_L4,absolute_difference")
    for row in rows
        @printf(
            "%s,%.16g,%.16g,%+.16g,%.16g\n",
            row.observable,
            row.L4,
            row.L100,
            row.signed_difference_L100_minus_L4,
            row.absolute_difference,
        )
    end
    @printf("validation_L4_identity_error=%.17g\n", L4.identity_error)
    @printf("validation_L100_identity_error=%.17g\n", L100.identity_error)
    @printf("validation_L4_density_reference_error=%.17g\n", L4_reference_errors.density)
    @printf("validation_L4_kinetic_reference_error=%.17g\n", L4_reference_errors.kinetic_energy_per_cell)
    @printf("validation_L4_grand_energy_reference_error=%.17g\n", L4_reference_errors.grand_energy_per_cell)
    @printf("validation_high_symmetry_max_abs_error=%.17g\n", high_symmetry_error)
    @printf("validation_band_path_high_symmetry_max_abs_error=%.17g\n", band_path_error)
    @printf("validation_L4_formula_max_abs_error=%.17g\n", maximum(Tuple(L4.formula_errors)))
    @printf("validation_L100_formula_max_abs_error=%.17g\n", maximum(Tuple(L100.formula_errors)))
end

function main(arguments=ARGS)
    (PARAMETERS.t1, PARAMETERS.t2, PARAMETERS.t3, PARAMETERS.t4) ==
        (-1.0, -1.75, -0.85, -0.65) || error("unexpected model parameters")
    generation_mode = prepare_output_directory(arguments)

    L4 = check_thermodynamics(L_SMALL)
    L100 = check_thermodynamics(L_LARGE)
    L4_reference_errors = check_l4_reference(L4.thermodynamics)
    high_symmetry_rows, high_symmetry_error = check_high_symmetry_points()
    band = band_path(
        points_per_segment=BAND_POINTS_PER_SEGMENT,
        parameters=PARAMETERS,
    )
    band_path_error = check_band_path_high_symmetry(band)
    rows = comparison_rows(L4.thermodynamics, L100.thermodynamics)
    all(row -> all(isfinite, (
        row.L4,
        row.L100,
        row.signed_difference_L100_minus_L4,
        row.absolute_difference,
    )), rows) || error("finite-size comparison contains a nonfinite value")

    momenta, lower_mesh, upper_mesh = visualization_mesh()
    configure_plot_defaults()
    band_plot = make_band_plot(band)
    fermi_surface_plot = make_fermi_surface_plot(
        momenta,
        lower_mesh,
        upper_mesh,
    )
    footer_plot = make_footer_plot()
    combined_layout = @layout [
        band_panel{0.58w} fermi_surface_panel
        footer_panel{0.08h}
    ]
    combined_plot = plot(
        band_plot,
        fermi_surface_plot,
        footer_plot;
        layout=combined_layout,
        size=(1800, 900),
        dpi=180,
        plot_title="Noninteracting two-orbital band structure and Fermi surface",
        plot_titlefontsize=16,
    )

    mkpath(OUTPUT_DIRECTORY)
    write_comparison(COMPARISON_PATH, rows)
    write_band_path(BAND_PATH_DATA_PATH, band)
    savefig(plot(band_plot; size=(1250, 820), dpi=180), BAND_FIGURE_PATH)
    savefig(
        plot(fermi_surface_plot; size=(1000, 900), dpi=180),
        FERMI_SURFACE_FIGURE_PATH,
    )
    savefig(combined_plot, COMBINED_FIGURE_PATH)
    write_metadata(
        METADATA_PATH,
        L4,
        L100,
        L4_reference_errors,
        high_symmetry_rows,
        high_symmetry_error,
        band_path_error,
        generation_mode,
    )

    for path in (
        BAND_FIGURE_PATH,
        FERMI_SURFACE_FIGURE_PATH,
        COMBINED_FIGURE_PATH,
        COMPARISON_PATH,
        BAND_PATH_DATA_PATH,
        METADATA_PATH,
    )
        isfile(path) && filesize(path) > 0 || error("missing output: $path")
        println("created=$path")
    end
    print_report(
        rows,
        L4,
        L100,
        L4_reference_errors,
        high_symmetry_error,
        band_path_error,
    )
    println("dqmc_state_constructed=false")
    println("metropolis_sampling_performed=false")
    return nothing
end

main()
