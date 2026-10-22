#!/usr/bin/env julia

# Deterministically extend the accepted 20-row FSS input with the validated
# L20 and L22 DQMC snapshots. The accepted input is read but never modified.

using SHA

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const BASE_INPUT = joinpath(REPO_ROOT, "results", "processed", "oneD_hubbard_fss_fit_input.csv")
const BASE_INPUT_SHA256 = "713cb6de4f5be017debf3d3a39701b47307d333d3908f2ca6106d914b74d1427"
const SNAPSHOT_COMMIT = "c766667ee7ee327bbf0744bfacbe01b82fd4cea0"

function usage(io::IO = stdout)
    println(io, "Usage:")
    println(io, "  julia analysis/build_oneD_hubbard_fss_L4_to_L22_input.jl --output OUTPUT.csv --no-clobber")
end

function parse_args(args::Vector{String})
    args == ["--help"] && (usage(); return nothing)
    length(args) == 3 && args[1] == "--output" && args[3] == "--no-clobber" ||
        error("Expected --output OUTPUT.csv --no-clobber; use --help")
    isempty(args[2]) && error("Output path must not be empty")
    return abspath(args[2])
end

sha256_file(path::String) = bytes2hex(sha256(read(path)))

function extended_row(header::Vector{String}, template::Vector{String};
        id::String, L::String, nA::String, beta::String,
        measurement_path::String, measurement_sha256::String, measurement_oid::String,
        parameter_path::String, parameter_sha256::String, parameter_oid::String,
        S2::String, sigma::String, raw_manifest_sha256::String)
    fields = copy(template)
    column = Dict(name => i for (i, name) in enumerate(header))
    values = Dict(
        "observation_id" => id,
        "measurement_source_path" => measurement_path,
        "measurement_source_sha256" => measurement_sha256,
        "measurement_source_git_blob_oid" => measurement_oid,
        "parameter_source_path" => parameter_path,
        "parameter_source_sha256" => parameter_sha256,
        "parameter_source_git_blob_oid" => parameter_oid,
        "source_commit" => SNAPSHOT_COMMIT,
        "L" => L,
        "nA" => nA,
        "beta_text" => beta,
        "S2_text" => S2,
        "S2_error_text" => sigma,
        "convergence_status" => "accepted_production_complete_and_dmrg_consistent",
        "notes" => join([
            "measurement_and_sigma_tokens_preserved_from_processed_summary",
            "production_commit=c4c873ffa16908a6641f851cb6c9ad41a9a22c29",
            "N_burnin=500",
            "N_measurements=4000",
            "N_bins=40",
            "raw_node_pairs=68",
            "raw_bin_rows=2720",
            "raw_combined_manifest_sha256=$raw_manifest_sha256",
        ], ";"),
    )
    for (name, value) in values
        fields[column[name]] = value
    end
    any(value -> occursin(',', value) || occursin('"', value) || occursin('\n', value) || occursin('\r', value), fields) &&
        error("Generated row unexpectedly requires CSV quoting")
    return join(fields, ',')
end

function main(args::Vector{String})
    output = parse_args(args)
    output === nothing && return 0
    isfile(BASE_INPUT) || error("Accepted base input is missing: $BASE_INPUT")
    sha256_file(BASE_INPUT) == BASE_INPUT_SHA256 || error("Accepted base-input SHA-256 mismatch")
    ispath(output) && error("Refusing to overwrite existing output: $output")

    lines = readlines(BASE_INPUT)
    length(lines) == 21 || error("Accepted base input must contain one header and 20 rows")
    any(line -> occursin('"', line), lines) && error("Accepted base input unexpectedly contains CSV quotes")
    header = String.(split(lines[1], ','; keepempty = true))
    insertion_index = findfirst(line -> startswith(line, "dqmc_L18,"), lines)
    insertion_index === nothing && error("Accepted dqmc_L18 row is missing")
    any(line -> startswith(line, "dqmc_L20,") || startswith(line, "dqmc_L22,"), lines) &&
        error("Accepted base input already contains L20 or L22")
    template = String.(split(lines[insertion_index], ','; keepempty = true))
    length(template) == length(header) || error("Template row width mismatch")

    l20 = extended_row(header, template;
        id = "dqmc_L20", L = "20", nA = "10", beta = "80",
        measurement_path = "results/processed/FSS_1D_U2_beta4L_L20_Ninc16.csv",
        measurement_sha256 = "d6c734257a969f6ccc465972fc7c079fa2772b0c9d9794820538a902db6b048b",
        measurement_oid = "47cec81f8eae2d5196113b23f54b270671383d7f",
        parameter_path = "parameters/L20_U2_beta80_nA10_Ninc16_highstat.toml",
        parameter_sha256 = "3096a974aa7f1325d9ee00b00933964112217605135a26ecdd5cdebcbcdd6803",
        parameter_oid = "32814c8d6270dcf571ab46c68e62a4d63f0885e6",
        S2 = "1.805179361272e+00", sigma = "8.209270540022e-03",
        raw_manifest_sha256 = "1d9e9fd097352a9cccf404a43d973abdca1aadc40d618dd5727c9024708d1609")
    l22 = extended_row(header, template;
        id = "dqmc_L22", L = "22", nA = "11", beta = "88",
        measurement_path = "results/processed/FSS_1D_U2_beta4L_L22_Ninc16.csv",
        measurement_sha256 = "bf82b619ae420912ec3795731b34c2791f3e4c37258a9195cadb8179cbe3dafd",
        measurement_oid = "9cfd7a90cb7962e970937e0f0019f3a38dc7e14b",
        parameter_path = "parameters/L22_U2_beta88_nA11_Ninc16_highstat.toml",
        parameter_sha256 = "e258e16c021eddbe0ede9473fcfd0cd0aba7eed408d460190b7c8203641dab75",
        parameter_oid = "bce7e2213a9698da97d5f01511c33e5071806ac9",
        S2 = "2.025158749526e+00", sigma = "4.730458606235e-03",
        raw_manifest_sha256 = "368ec6bcce670690c3105192c4afb58926a4c38a377378611335b83926c09eca")

    output_lines = vcat(lines[1:insertion_index], [l20, l22], lines[insertion_index + 1:end])
    mkpath(dirname(output))
    write(output, join(output_lines, '\n') * "\n")
    println("Wrote 22 canonical observations: $output")
    println("SHA-256: ", sha256_file(output))
    return 0
end

try
    exit(main(ARGS))
catch error_object
    println(stderr, "ERROR: ", sprint(showerror, error_object))
    exit(2)
end
