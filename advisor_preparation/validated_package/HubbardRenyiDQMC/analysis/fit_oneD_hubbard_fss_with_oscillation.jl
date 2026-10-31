#!/usr/bin/env julia

# Reproducible finite-size diagnostic fits for the one-dimensional U/t=2 data.
# The reported c_eff values are finite-size regression diagnostics only; they are
# not asymptotic central-charge determinations or physical uncertainty estimates.

ENV["GKSwstype"] = "100"

using LinearAlgebra
using Printf
using Plots

const LEGACY_SOURCE_COMMIT = "88c21ae3345ad7414d86b37e777e10c35e6451f2"
const EXTENDED_DQMC_SOURCE_COMMIT = "c766667ee7ee327bbf0744bfacbe01b82fd4cea0"
const EXTENDED_DQMC_IDS = Set(["dqmc_L20", "dqmc_L22"])
const REQUIRED_COLUMNS = [
    "observation_id",
    "method",
    "measurement_source_path",
    "measurement_source_selector",
    "measurement_source_sha256",
    "measurement_source_git_blob_oid",
    "parameter_source_path",
    "parameter_source_sha256",
    "parameter_source_git_blob_oid",
    "source_commit",
    "L",
    "t_text",
    "U_text",
    "U_over_t_text",
    "mu_text",
    "ensemble",
    "filling",
    "Nup_text",
    "Ndn_text",
    "boundary_condition",
    "nA",
    "nA_convention",
    "subsystem_sites",
    "renyi_order",
    "renyi_definition",
    "state_convention",
    "beta_text",
    "beta_rule",
    "S2_text",
    "S2_error_text",
    "uncertainty_type",
    "S2_source_significant_digits",
    "S2_error_source_significant_digits",
    "convergence_status",
    "overlap_status",
    "seed_text",
    "sweeps_text",
    "maxdim_text",
    "maxlinkdim_text",
    "cutoff_text",
    "truncation_error_text",
    "energy_density_text",
    "normcheck_text",
    "notes",
]

const EXPECTED_DQMC_L_SETS = (
    [4, 6, 8, 10, 12, 14, 16, 18],
    [4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
)

const EXPECTED_REFERENCE_L = Dict(
    "ED" => [4, 6, 8, 10, 12],
    "DMRG" => [14, 16, 18, 20, 22, 24, 26],
)

const BRANCH_0 = "L mod 4 = 0"
const BRANCH_2 = "L mod 4 = 2"
const ALL_BRANCH = "all"
const DATASET_ORDER = ["dqmc", "reference_without_L26", "reference_with_L26_candidate"]
const ABS_TOL = 1.0e-10
const REL_TOL = 1.0e-9
const MAX_NORMALIZED_CONDITION = 1.0e6
const WARN_NORMALIZED_CONDITION = 1.0e4

struct Observation
    id::String
    method::String
    L::Int
    nA::Int
    S2::Float64
    sigma::Union{Nothing,Float64}
    uncertainty_type::String
    beta_text::String
    source_row::Dict{String,String}
end

struct FitResult
    fit_id::String
    dataset_variant::String
    role::String
    model::String
    branch::String
    p::Union{Nothing,Int}
    L_min::Int
    L_max::Int
    window_kind::String
    omitted_observation_id::String
    n_points::Int
    n_parameters::Int
    dof::Int
    rank::Int
    condition_number_raw::Float64
    condition_number_normalized::Float64
    condition_warning::Bool
    c_eff_diagnostic::Float64
    c_eff_nominal_se::Float64
    b::Float64
    b_nominal_se::Float64
    A::Union{Nothing,Float64}
    A_nominal_se::Union{Nothing,Float64}
    se_kind::String
    statistic_name::String
    statistic_value::Float64
    input_sha256::String
    l26_status::String
    qr_svd_max_abs_delta::Float64
    notes::String
    coefficients::Vector{Float64}
end

struct Options
    input::String
    summary_output::String
    figure_dir::String
    no_clobber::Bool
end

function usage(io::IO = stdout)
    println(io, "Usage:")
    println(io, "  julia --project=ENV analysis/fit_oneD_hubbard_fss_with_oscillation.jl \\")
    println(io, "    --input INPUT.csv --summary-output SUMMARY.csv \\")
    println(io, "    --figure-dir FIGURE_DIR --no-clobber")
    println(io)
    println(io, "Required options:")
    println(io, "  --input PATH           canonical 20- or 22-row provenance CSV")
    println(io, "  --summary-output PATH  new fit-summary CSV")
    println(io, "  --figure-dir PATH      directory for four new PNG/SVG pairs")
    println(io, "  --no-clobber           refuse every existing output target")
    println(io, "  --help                 show this help and exit")
end

function parse_options(args::Vector{String})
    isempty(args) && error("Missing required arguments; use --help")
    args == ["--help"] && (usage(); return nothing)

    values = Dict{String,String}()
    no_clobber = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help"
            error("--help must be used by itself")
        elseif arg == "--no-clobber"
            no_clobber && error("Duplicate option: --no-clobber")
            no_clobber = true
            i += 1
        elseif arg in ("--input", "--summary-output", "--figure-dir")
            haskey(values, arg) && error("Duplicate option: $arg")
            i == length(args) && error("Missing value for $arg")
            startswith(args[i + 1], "--") && error("Missing value for $arg")
            isempty(args[i + 1]) && error("Empty value for $arg")
            values[arg] = args[i + 1]
            i += 2
        else
            error("Unknown argument: $arg")
        end
    end

    for arg in ("--input", "--summary-output", "--figure-dir")
        haskey(values, arg) || error("Missing required option: $arg")
    end
    no_clobber || error("--no-clobber is required for non-destructive generation")
    return Options(abspath(values["--input"]), abspath(values["--summary-output"]), abspath(values["--figure-dir"]), true)
end

# Minimal RFC 4180 parser: quotes, doubled quotes, CRLF, and commas are supported.
# Embedded newlines are supported too, although canonical inputs do not use them.
function parse_csv(text::String)
    chars = collect(text)
    rows = Vector{Vector{String}}()
    row = String[]
    field = IOBuffer()
    quoted = false
    i = 1
    while i <= length(chars)
        c = chars[i]
        if quoted
            if c == '"'
                if i < length(chars) && chars[i + 1] == '"'
                    write(field, '"')
                    i += 2
                    continue
                end
                quoted = false
            else
                write(field, c)
            end
        else
            if c == '"'
                position(field) == 0 || error("Malformed CSV: quote inside unquoted field")
                quoted = true
            elseif c == ','
                push!(row, String(take!(field)))
            elseif c == '\n'
                push!(row, String(take!(field)))
                push!(rows, row)
                row = String[]
            elseif c == '\r'
                if i == length(chars) || chars[i + 1] != '\n'
                    push!(row, String(take!(field)))
                    push!(rows, row)
                    row = String[]
                end
            else
                write(field, c)
            end
        end
        i += 1
    end
    quoted && error("Malformed CSV: unterminated quoted field")
    if position(field) > 0 || !isempty(row)
        push!(row, String(take!(field)))
        push!(rows, row)
    end
    while !isempty(rows) && all(isempty, rows[end])
        pop!(rows)
    end
    isempty(rows) && error("Empty CSV")
    return rows
end

function read_csv_dicts(path::String)
    isfile(path) || error("Input is not a file: $path")
    rows = parse_csv(read(path, String))
    header = rows[1]
    isempty(header) && error("CSV header is empty")
    length(unique(header)) == length(header) || error("CSV header contains duplicate columns")
    missing = setdiff(REQUIRED_COLUMNS, header)
    isempty(missing) || error("CSV is missing required columns: $(join(missing, ", "))")
    result = Dict{String,String}[]
    for (offset, fields) in enumerate(rows[2:end])
        line_number = offset + 1
        length(fields) == length(header) || error("CSV row $line_number has $(length(fields)) fields; expected $(length(header))")
        all(isempty, fields) && continue
        push!(result, Dict(header[i] => fields[i] for i in eachindex(header)))
    end
    return result
end

function parse_finite_float(text::String, label::String)
    isempty(strip(text)) && error("Missing numeric value for $label")
    value = tryparse(Float64, text)
    value === nothing && error("Malformed numeric value for $label: $text")
    isfinite(value) || error("Non-finite numeric value for $label: $text")
    return value
end

function parse_int_exact(text::String, label::String)
    isempty(strip(text)) && error("Missing integer value for $label")
    value = tryparse(Int, text)
    value === nothing && error("Malformed integer value for $label: $text")
    return value
end

function require_close(actual::Real, expected::Real, label::String)
    tolerance = ABS_TOL + REL_TOL * max(abs(actual), abs(expected))
    abs(actual - expected) <= tolerance || error("$label mismatch: expected $expected, got $actual (tolerance $tolerance)")
end

is_sha256(text::String) = occursin(r"^[0-9a-f]{64}$", text)
is_git_oid_or_na(text::String) = isempty(text) || text == "not_available_untracked_source" || occursin(r"^[0-9a-f]{40}$", text) || occursin(r"^[0-9a-f]{64}$", text)

branch_label(L::Int) = mod(L, 4) == 0 ? BRANCH_0 : BRANCH_2

function validate_input_rows(raw_rows::Vector{Dict{String,String}})
    length(raw_rows) in (20, 22) || error("Canonical input must contain exactly 20 or 22 data rows; found $(length(raw_rows))")
    ids = [row["observation_id"] for row in raw_rows]
    all(!isempty, ids) || error("Every row requires a nonempty observation_id")
    length(unique(ids)) == length(ids) || error("Duplicate observation_id in canonical input")

    for (offset, row) in enumerate(raw_rows)
        i = offset + 1
        id = row["observation_id"]
        method = row["method"]
        method == "DQMC" || haskey(EXPECTED_REFERENCE_L, method) || error("Row $i ($id): unexpected method $method")
        L = parse_int_exact(row["L"], "$id L")
        nA = parse_int_exact(row["nA"], "$id nA")
        iseven(L) || error("$id: L must be even")
        nA == div(L, 2) || error("$id: nA must equal L/2")
        parse_int_exact(row["renyi_order"], "$id renyi_order") == 2 || error("$id: renyi_order must be 2")
        parse_int_exact(row["S2_source_significant_digits"], "$id S2 source digits") > 0 || error("$id: S2 source significant digits must be positive")

        require_close(parse_finite_float(row["t_text"], "$id t"), 1.0, "$id t")
        require_close(parse_finite_float(row["U_text"], "$id U"), 2.0, "$id U")
        require_close(parse_finite_float(row["U_over_t_text"], "$id U/t"), 2.0, "$id U/t")
        parse_finite_float(row["S2_text"], "$id S2")

        expected_source_commit = method == "DQMC" && L in (20, 22) ? EXTENDED_DQMC_SOURCE_COMMIT : LEGACY_SOURCE_COMMIT
        row["source_commit"] == expected_source_commit || error("$id: source_commit mismatch")
        isempty(row["measurement_source_path"]) && error("$id: missing measurement source path")
        isempty(row["measurement_source_selector"]) && error("$id: missing measurement source selector")
        is_sha256(row["measurement_source_sha256"]) || error("$id: invalid measurement source SHA-256")
        is_git_oid_or_na(row["measurement_source_git_blob_oid"]) || error("$id: invalid measurement source git blob/OID status")
        isempty(row["parameter_source_path"]) && error("$id: missing parameter source path")
        is_sha256(row["parameter_source_sha256"]) || error("$id: invalid parameter source SHA-256")
        is_git_oid_or_na(row["parameter_source_git_blob_oid"]) || error("$id: invalid parameter source git blob/OID status")

        occursin("half", lowercase(row["filling"])) || error("$id: filling must explicitly identify half filling")
        startswith(lowercase(row["boundary_condition"]), "periodic") || error("$id: boundary condition must be periodic")
        isempty(row["nA_convention"]) && error("$id: nA convention is missing")
        row["subsystem_sites"] in ("1:nA", "1:$nA") || error("$id: subsystem_sites must encode sites 1:nA")
        isempty(row["renyi_definition"]) && error("$id: Renyi definition is missing")
        isempty(row["ensemble"]) && error("$id: ensemble is missing")
        isempty(row["state_convention"]) && error("$id: state convention is missing")
        isempty(row["beta_rule"]) && error("$id: beta rule is missing")
        isempty(row["uncertainty_type"]) && error("$id: uncertainty type is missing")
        isempty(row["convergence_status"]) && error("$id: convergence status is missing")

        if method == "DQMC"
            lowercase(row["ensemble"]) == "grand_canonical" || error("$id: DQMC ensemble must be grand_canonical")
            require_close(parse_finite_float(row["mu_text"], "$id mu"), 0.0, "$id mu")
            occursin("finite", lowercase(row["state_convention"])) || error("$id: DQMC state must be finite beta")
            require_close(parse_finite_float(row["beta_text"], "$id beta"), 4.0L, "$id beta=4L")
            occursin("4", row["beta_rule"]) && occursin("L", row["beta_rule"]) || error("$id: DQMC beta_rule must encode beta=4L")
            sigma = parse_finite_float(row["S2_error_text"], "$id S2 error")
            sigma > 0.0 || error("$id: DQMC S2 error must be positive")
            row["uncertainty_type"] == "independent_absolute_mc_1sigma" || error("$id: DQMC uncertainty_type must be independent_absolute_mc_1sigma")
            parse_int_exact(row["S2_error_source_significant_digits"], "$id S2 error digits") > 0 || error("$id: S2 error digits must be positive")
            isempty(row["Nup_text"]) || error("$id: grand-canonical DQMC Nup_text must be empty")
            isempty(row["Ndn_text"]) || error("$id: grand-canonical DQMC Ndn_text must be empty")
        else
            lowercase(row["ensemble"]) in ("fixed_particle_number", "canonical_fixed_particle_number") || error("$id: reference ensemble must be canonical fixed particle number")
            occursin("ground", lowercase(row["state_convention"])) || error("$id: reference state must be ground state")
            lowercase(row["beta_text"]) in ("inf", "infinity", "ground_state") || error("$id: reference beta_text must identify the ground-state limit")
            isempty(row["S2_error_text"]) || error("$id: deterministic reference must not fabricate S2_error")
            isempty(row["S2_error_source_significant_digits"]) || error("$id: deterministic reference must leave error significant digits empty")
            parse_int_exact(row["Nup_text"], "$id Nup") == div(L, 2) || error("$id: Nup must equal L/2")
            parse_int_exact(row["Ndn_text"], "$id Ndn") == div(L, 2) || error("$id: Ndn must equal L/2")
            if method == "ED"
                row["uncertainty_type"] == "deterministic_exact_no_statistical_error" || error("$id: unexpected ED uncertainty type")
            else
                row["uncertainty_type"] == "not_reported_convergence_limited" || error("$id: unexpected DMRG uncertainty type")
            end
        end
    end

    observed_dqmc = sort([parse_int_exact(row["L"], "L") for row in raw_rows if row["method"] == "DQMC"])
    observed_dqmc in EXPECTED_DQMC_L_SETS || error("DQMC L set mismatch: got $observed_dqmc")
    extended_ids_observed = Set(
        row["observation_id"] for row in raw_rows
        if row["method"] == "DQMC" && row["L"] in ("20", "22")
    )
    if observed_dqmc == EXPECTED_DQMC_L_SETS[2]
        extended_ids_observed == EXTENDED_DQMC_IDS || error("Extended DQMC observation IDs must be dqmc_L20 and dqmc_L22")
    else
        isempty(intersect(Set(ids), EXTENDED_DQMC_IDS)) || error("Legacy DQMC input cannot reuse the L20/L22 observation IDs")
    end
    for method in ("ED", "DMRG")
        observed = sort([parse_int_exact(row["L"], "L") for row in raw_rows if row["method"] == method])
        observed == EXPECTED_REFERENCE_L[method] || error("$method L set mismatch: expected $(EXPECTED_REFERENCE_L[method]), got $observed")
    end

    conventions = unique([(row["nA_convention"], row["renyi_definition"]) for row in raw_rows])
    length(conventions) == 1 || error("All rows must use one nA convention and one Renyi definition")

    ed10 = only(row for row in raw_rows if row["method"] == "ED" && row["L"] == "10")
    ed10["S2_text"] == "1.716551661433" || error("ED L10 must preserve repository-authoritative token 1.716551661433")

    l26 = only(row for row in raw_rows if row["method"] == "DMRG" && row["L"] == "26")
    l26["S2_text"] == "2.0852574848816734" || error("L26 canonical S2 token mismatch")
    l26["energy_density_text"] == "-0.84679839487120911" || error("L26 energy-density token mismatch")
    l26["seed_text"] == "1234" || error("L26 canonical seed must be 1234")
    l26["sweeps_text"] == "22" || error("L26 canonical sweep count must be 22")
    l26["maxdim_text"] == "5600" || error("L26 canonical maxdim must be 5600")
    l26["maxlinkdim_text"] == "5600" || error("L26 canonical maxlinkdim must be 5600")
    l26["overlap_status"] == "not_available_at_HEAD" || error("L26 overlap status mismatch")
    truncation = parse_finite_float(l26["truncation_error_text"], "L26 truncation error")
    truncation <= 1.0e-8 || error("L26 truncation error exceeds 1e-8")
    normcheck = parse_finite_float(l26["normcheck_text"], "L26 normcheck")
    abs(1.0 - normcheck) <= 1.0e-10 || error("L26 normcheck fails 1e-10 gate")

    return nothing
end

function to_observations(raw_rows::Vector{Dict{String,String}})
    observations = Observation[]
    for row in raw_rows
        method = row["method"]
        sigma = method == "DQMC" ? parse(Float64, row["S2_error_text"]) : nothing
        push!(observations, Observation(
            row["observation_id"], method, parse(Int, row["L"]), parse(Int, row["nA"]),
            parse(Float64, row["S2_text"]), sigma, row["uncertainty_type"], row["beta_text"], row,
        ))
    end
    sort!(observations; by = observation -> (findfirst(==(observation.method), ["DQMC", "ED", "DMRG"]), observation.L))
    return observations
end

# Pure-Julia SHA-256 keeps the analysis environment limited to the declared
# LinearAlgebra, Printf, and Plots dependencies.
const SHA256_K = UInt32[
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]

rotr(x::UInt32, n::Int) = (x >> n) | (x << (32 - n))

function sha256_hex(bytes::Vector{UInt8})
    message = copy(bytes)
    bit_length = UInt64(length(message)) * UInt64(8)
    push!(message, 0x80)
    while mod(length(message), 64) != 56
        push!(message, 0x00)
    end
    for shift in 56:-8:0
        push!(message, UInt8((bit_length >> shift) & 0xff))
    end

    h = UInt32[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
    w = Vector{UInt32}(undef, 64)
    for offset in 1:64:length(message)
        for i in 0:15
            j = offset + 4i
            w[i + 1] = (UInt32(message[j]) << 24) | (UInt32(message[j + 1]) << 16) | (UInt32(message[j + 2]) << 8) | UInt32(message[j + 3])
        end
        for i in 17:64
            s0 = rotr(w[i - 15], 7) ⊻ rotr(w[i - 15], 18) ⊻ (w[i - 15] >> 3)
            s1 = rotr(w[i - 2], 17) ⊻ rotr(w[i - 2], 19) ⊻ (w[i - 2] >> 10)
            w[i] = w[i - 16] + s0 + w[i - 7] + s1
        end
        a,b,c,d,e,f,g,hh = h
        for i in 1:64
            s1 = rotr(e, 6) ⊻ rotr(e, 11) ⊻ rotr(e, 25)
            ch = (e & f) ⊻ ((~e) & g)
            temp1 = hh + s1 + ch + SHA256_K[i] + w[i]
            s0 = rotr(a, 2) ⊻ rotr(a, 13) ⊻ rotr(a, 22)
            maj = (a & b) ⊻ (a & c) ⊻ (b & c)
            temp2 = s0 + maj
            hh,g,f,e,d,c,b,a = g,f,e,d + temp1,c,b,a,temp1 + temp2
        end
        h .= h .+ UInt32[a,b,c,d,e,f,g,hh]
    end
    return join(@sprintf("%08x", value) for value in h)
end

sha256_file(path::String) = sha256_hex(read(path))

function design_matrix(points::Vector{Observation}, model::String, p::Union{Nothing,Int})
    x = [0.25 * log(point.L / pi) for point in points]
    if model in ("total", "branch")
        return hcat(x, ones(length(points)))
    elseif model == "oscillation"
        p in (1, 2) || error("Oscillation model requires p=1 or p=2")
        oscillator = [(-1.0)^div(point.L, 2) / point.L^p for point in points]
        return hcat(x, ones(length(points)), oscillator)
    end
    error("Unknown model: $model")
end

function svd_rank_and_condition(matrix::Matrix{Float64})
    decomposition = svd(matrix)
    singular_values = decomposition.S
    threshold = maximum(size(matrix)) * eps(Float64) * maximum(singular_values)
    rank_value = count(>(threshold), singular_values)
    condition = singular_values[1] / singular_values[end]
    return decomposition, rank_value, condition
end

function fit_model(
    dataset_variant::String,
    role::String,
    model::String,
    branch::String,
    p::Union{Nothing,Int},
    selected::Vector{Observation},
    window_kind::String,
    omitted_id::String,
    input_sha256::String,
    sequence_number::Int,
)
    points = sort(selected; by = point -> point.L)
    X = Matrix{Float64}(design_matrix(points, model, p))
    y = Float64[point.S2 for point in points]
    n_points, n_parameters = size(X)
    n_points > n_parameters || error("Underdetermined fit: $dataset_variant/$role/$model/$branch uses n=$n_points, p=$n_parameters")
    dof = n_points - n_parameters
    dof > 0 || error("Nonpositive degrees of freedom for $dataset_variant/$role/$model/$branch")

    weighted = dataset_variant == "dqmc"
    if weighted
        sigma = Float64[something(point.sigma) for point in points]
        all(isfinite, sigma) && all(>(0.0), sigma) || error("DQMC fit contains non-finite or nonpositive sigma")
        solve_matrix = X ./ sigma
        solve_y = y ./ sigma
    else
        solve_matrix = X
        solve_y = y
    end

    column_norms = [norm(view(solve_matrix, :, j)) for j in axes(solve_matrix, 2)]
    all(isfinite, column_norms) && all(>(0.0), column_norms) || error("Degenerate design-matrix column")
    normalized = solve_matrix ./ reshape(column_norms, 1, :)
    decomposition, rank_value, raw_condition = svd_rank_and_condition(solve_matrix)
    _, normalized_rank, normalized_condition = svd_rank_and_condition(normalized)
    rank_value == n_parameters || error("Rank-deficient fit: $dataset_variant/$role/$model/$branch rank=$rank_value, p=$n_parameters")
    normalized_rank == n_parameters || error("Column-normalized matrix is rank deficient")
    isfinite(normalized_condition) || error("Non-finite normalized condition number")
    normalized_condition < MAX_NORMALIZED_CONDITION || error("Normalized condition number $normalized_condition exceeds strict limit $MAX_NORMALIZED_CONDITION")

    coefficients_qr = qr(solve_matrix) \ solve_y
    coefficients_svd = decomposition.V * ((transpose(decomposition.U) * solve_y) ./ decomposition.S)
    all(isfinite, coefficients_qr) || error("QR returned non-finite coefficients")
    all(isfinite, coefficients_svd) || error("SVD returned non-finite coefficients")
    deltas = abs.(coefficients_qr .- coefficients_svd)
    for j in eachindex(coefficients_qr)
        tolerance = ABS_TOL + REL_TOL * max(abs(coefficients_qr[j]), abs(coefficients_svd[j]))
        deltas[j] <= tolerance || error("QR/SVD coefficient mismatch at column $j: $(coefficients_qr[j]) vs $(coefficients_svd[j])")
    end

    residual = y - X * coefficients_qr
    inverse_normal = decomposition.V * Diagonal(1.0 ./ (decomposition.S .^ 2)) * transpose(decomposition.V)
    if weighted
        statistic_name = "chi2_per_dof"
        statistic_value = sum(abs2, residual ./ Float64[something(point.sigma) for point in points]) / dof
        covariance = inverse_normal
        se_kind = "nominal_wls_se"
        fit_note = "DQMC errors treated as independent absolute 1sigma; Sigma=diag(sigma_i^2); covariance=(X'WX)^-1 without reduced-chi2 rescaling"
    else
        statistic_name = "rss_per_dof"
        statistic_value = sum(abs2, residual) / dof
        covariance = statistic_value .* inverse_normal
        se_kind = "nominal_ols_se"
        fit_note = "nominal OLS regression scatter only; not ED precision, DMRG uncertainty, or physical uncertainty"
    end
    standard_errors = sqrt.(max.(0.0, diag(covariance)))
    all(isfinite, standard_errors) || error("Non-finite nominal regression standard error")

    l26_status = dataset_variant == "reference_with_L26_candidate" ? "provisional_candidate_included" :
                  dataset_variant == "reference_without_L26" ? "excluded" : "not_applicable"
    ptext = p === nothing ? "na" : string(p)
    fit_id = @sprintf("fit_%03d_%s_%s_%s_%s_p%s", sequence_number, dataset_variant, role, model, replace(branch, " " => "_", "=" => "eq"), ptext)
    condition_note = normalized_condition >= WARN_NORMALIZED_CONDITION ? "; normalized condition warning >=1e4" : ""
    A = n_parameters == 3 ? coefficients_qr[3] : nothing
    A_se = n_parameters == 3 ? standard_errors[3] : nothing
    return FitResult(
        fit_id, dataset_variant, role, model, branch, p,
        minimum(point.L for point in points), maximum(point.L for point in points),
        window_kind, omitted_id, n_points, n_parameters, dof, rank_value,
        raw_condition, normalized_condition, normalized_condition >= WARN_NORMALIZED_CONDITION,
        coefficients_qr[1], standard_errors[1], coefficients_qr[2], standard_errors[2],
        A, A_se, se_kind, statistic_name, statistic_value, input_sha256, l26_status,
        maximum(deltas),
        "c_eff is a finite-size diagnostic only; $fit_note$condition_note",
        collect(coefficients_qr),
    )
end

function dataset_points(observations::Vector{Observation}, dataset_variant::String)
    if dataset_variant == "dqmc"
        return sort([point for point in observations if point.method == "DQMC"]; by = point -> point.L)
    elseif dataset_variant == "reference_without_L26"
        return sort([point for point in observations if point.method in ("ED", "DMRG") && point.L != 26]; by = point -> point.L)
    elseif dataset_variant == "reference_with_L26_candidate"
        return sort([point for point in observations if point.method in ("ED", "DMRG")]; by = point -> point.L)
    end
    error("Unknown dataset variant: $dataset_variant")
end

function run_all_fits(observations::Vector{Observation}, input_sha256::String)
    fits = FitResult[]
    sequence = Ref(0)
    function add_fit(dataset, role, model, branch, p, points, window_kind; omitted = "")
        sequence[] += 1
        push!(fits, fit_model(dataset, role, model, branch, p, points, window_kind, omitted, input_sha256, sequence[]))
    end

    # Principal full-window diagnostics.
    for dataset in DATASET_ORDER
        points = dataset_points(observations, dataset)
        add_fit(dataset, "principal", "total", ALL_BRANCH, nothing, points, "full")
        for branch in (BRANCH_0, BRANCH_2)
            selected = [point for point in points if branch_label(point.L) == branch]
            add_fit(dataset, "principal", "branch", branch, nothing, selected, "full")
        end
        for p in (1, 2)
            add_fit(dataset, "principal", "oscillation", ALL_BRANCH, p, points, "full")
        end
    end

    # Approved Lmin sensitivity checks.
    for dataset in DATASET_ORDER
        points = dataset_points(observations, dataset)
        all_lmins = dataset == "dqmc" ? [4, 6, 8, 10] : [4, 6, 8, 10, 12]
        for Lmin in all_lmins
            selected = [point for point in points if point.L >= Lmin]
            add_fit(dataset, "lmin_sensitivity", "total", ALL_BRANCH, nothing, selected, "Lmin")
            for p in (1, 2)
                add_fit(dataset, "lmin_sensitivity", "oscillation", ALL_BRANCH, p, selected, "Lmin")
            end
        end
        branch_lmins = dataset == "dqmc" ? Dict(BRANCH_0 => [4, 8], BRANCH_2 => [6, 10]) :
                                            Dict(BRANCH_0 => [4, 8, 12], BRANCH_2 => [6, 10, 14])
        for branch in (BRANCH_0, BRANCH_2), Lmin in branch_lmins[branch]
            selected = [point for point in points if branch_label(point.L) == branch && point.L >= Lmin]
            add_fit(dataset, "lmin_sensitivity", "branch", branch, nothing, selected, "Lmin")
        end
    end

    # Every contiguous five-size all-branch window and three-point branch window.
    for dataset in DATASET_ORDER
        points = dataset_points(observations, dataset)
        for start in 1:(length(points) - 4)
            selected = points[start:start + 4]
            add_fit(dataset, "contiguous_window", "total", ALL_BRANCH, nothing, selected, "contiguous_5_size")
            for p in (1, 2)
                add_fit(dataset, "contiguous_window", "oscillation", ALL_BRANCH, p, selected, "contiguous_5_size")
            end
        end
        for branch in (BRANCH_0, BRANCH_2)
            branch_points = [point for point in points if branch_label(point.L) == branch]
            for start in 1:(length(branch_points) - 2)
                selected = branch_points[start:start + 2]
                add_fit(dataset, "contiguous_window", "branch", branch, nothing, selected, "contiguous_3_branch")
            end
        end
    end

    # Leave one observation out from every full-window baseline model. No LOO
    # is applied to minimal windows, so omission cannot underdetermine a fit.
    for dataset in DATASET_ORDER
        points = dataset_points(observations, dataset)
        for model in ("total", "oscillation"), p in (model == "total" ? (nothing,) : (1, 2))
            for omitted in points
                selected = [point for point in points if point.id != omitted.id]
                add_fit(dataset, "leave_one_out", model, ALL_BRANCH, p, selected, "full_leave_one_out"; omitted = omitted.id)
            end
        end
        for branch in (BRANCH_0, BRANCH_2)
            branch_points = [point for point in points if branch_label(point.L) == branch]
            for omitted in branch_points
                selected = [point for point in branch_points if point.id != omitted.id]
                add_fit(dataset, "leave_one_out", "branch", branch, nothing, selected, "full_leave_one_out"; omitted = omitted.id)
            end
        end
    end

    dqmc_count = count(point -> point.method == "DQMC", observations)
    expected_fit_count = dqmc_count == 8 ? 273 : dqmc_count == 10 ? 289 : error("Unexpected DQMC observation count: $dqmc_count")
    length(fits) == expected_fit_count || error("Internal fit-coverage mismatch: expected $expected_fit_count fits, generated $(length(fits))")
    return fits
end

function csv_field(value)
    text = if value === nothing
        ""
    elseif value isa Bool
        value ? "true" : "false"
    elseif value isa AbstractFloat
        isfinite(value) || error("Refusing to serialize non-finite value")
        @sprintf("%.17g", value)
    else
        string(value)
    end
    if occursin(',', text) || occursin('"', text) || occursin('\n', text) || occursin('\r', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

const SUMMARY_HEADER = [
    "fit_id", "dataset_variant", "role", "model", "branch", "p", "L_min", "L_max",
    "window_kind", "omitted_observation_id", "n_points", "n_parameters", "dof", "rank",
    "condition_number_raw", "condition_number_normalized", "condition_warning",
    "c_eff_diagnostic", "c_eff_nominal_se", "b", "b_nominal_se", "A", "A_nominal_se",
    "se_kind", "statistic_name", "statistic_value", "input_sha256", "l26_status",
    "qr_svd_max_abs_delta", "notes",
]

function summary_row(fit::FitResult)
    return Any[
        fit.fit_id, fit.dataset_variant, fit.role, fit.model, fit.branch, fit.p,
        fit.L_min, fit.L_max, fit.window_kind, fit.omitted_observation_id,
        fit.n_points, fit.n_parameters, fit.dof, fit.rank,
        fit.condition_number_raw, fit.condition_number_normalized, fit.condition_warning,
        fit.c_eff_diagnostic, fit.c_eff_nominal_se, fit.b, fit.b_nominal_se,
        fit.A, fit.A_nominal_se, fit.se_kind, fit.statistic_name, fit.statistic_value,
        fit.input_sha256, fit.l26_status, fit.qr_svd_max_abs_delta, fit.notes,
    ]
end

function write_summary(path::String, fits::Vector{FitResult})
    exclusive_file = Base.Filesystem.open(path, Base.JL_O_WRONLY | Base.JL_O_CREAT | Base.JL_O_EXCL, 0o644)
    close(exclusive_file)
    open(path, "w") do io
        println(io, join(SUMMARY_HEADER, ','))
        for fit in fits
            println(io, join(csv_field.(summary_row(fit)), ','))
        end
    end
end

function principal_fit(fits::Vector{FitResult}, dataset::String, model::String, branch::String; p = nothing)
    matches = [fit for fit in fits if fit.dataset_variant == dataset && fit.role == "principal" && fit.model == model && fit.branch == branch && fit.p == p]
    length(matches) == 1 || error("Expected one principal fit for $dataset/$model/$branch/p=$p")
    return only(matches)
end

function predict(fit::FitResult, L::Real)
    value = (fit.c_eff_diagnostic / 4.0) * log(L / pi) + fit.b
    if fit.model == "oscillation"
        integer_L = round(Int, L)
        value += something(fit.A) * ((-1.0)^div(integer_L, 2) / integer_L^something(fit.p))
    end
    return value
end

const COLOR_DQMC = "#1f77b4"
const COLOR_WITHOUT = "#2ca02c"
const COLOR_WITH = "#d62728"
const COLOR_ED = "#7f3c8d"
const COLOR_DMRG = "#11a579"
const COLOR_L26 = "#e68310"

function plotting_observations(observations::Vector{Observation}; branch::String = ALL_BRANCH)
    branch == ALL_BRANCH && return observations
    return [point for point in observations if branch_label(point.L) == branch]
end

function plot_points!(plt, observations::Vector{Observation}; branch::String = ALL_BRANCH, show_labels::Bool = false)
    points = plotting_observations(observations; branch = branch)
    dqmc = [point for point in points if point.method == "DQMC"]
    ed = [point for point in points if point.method == "ED"]
    dmrg = [point for point in points if point.method == "DMRG" && point.L != 26]
    l26 = [point for point in points if point.method == "DMRG" && point.L == 26]

    if !isempty(dqmc)
        scatter!(plt, [log(p.L / pi) for p in dqmc], [p.S2 for p in dqmc];
            yerror = [something(p.sigma) for p in dqmc], marker = :circle, markersize = 5,
            markercolor = COLOR_DQMC, linecolor = COLOR_DQMC,
            label = "DQMC finite beta=4L (reported 1sigma; diagonal fit)")
    end
    if !isempty(ed)
        scatter!(plt, [log(p.L / pi) for p in ed], [p.S2 for p in ed];
            marker = :xcross, markersize = 7, markerstrokewidth = 2,
            markercolor = COLOR_ED, label = "ED ground-state reference")
    end
    if !isempty(dmrg)
        scatter!(plt, [log(p.L / pi) for p in dmrg], [p.S2 for p in dmrg];
            marker = :square, markersize = 5, markercolor = COLOR_DMRG,
            label = "DMRG ground-state reference (L<=24)")
    end
    if !isempty(l26)
        scatter!(plt, [log(p.L / pi) for p in l26], [p.S2 for p in l26];
            marker = :star5, markersize = 10, markercolor = COLOR_L26,
            markerstrokecolor = :black, label = "L26 provisional DMRG candidate")
    end
    if show_labels
        for point in points
            annotate!(plt, log(point.L / pi), point.S2 + 0.025, text("L=$(point.L)", 7, :center))
        end
    end
    return plt
end

function base_plot(; title::String, size = (1200, 800), legend = :topleft)
    return plot(
        xlabel = "log(L/pi)", ylabel = "Renyi-2 entropy S2", title = title,
        legend = legend, size = size, dpi = 180, grid = true,
        framestyle = :box, foreground_color_legend = :black,
        background_color_legend = :white, legendfontsize = 8,
        guidefontsize = 11, tickfontsize = 9, titlefontsize = 12,
        left_margin = 8Plots.mm, right_margin = 6Plots.mm,
        top_margin = 6Plots.mm, bottom_margin = 7Plots.mm,
    )
end

function add_fit_line!(plt, fit::FitResult; label::String, color, linestyle = :solid, discrete::Bool = false)
    if discrete
        Ls = collect(fit.L_min:2:fit.L_max)
        xs = log.(Ls ./ pi)
        ys = [predict(fit, L) for L in Ls]
    else
        xs = range(log(fit.L_min / pi), log(fit.L_max / pi); length = 240)
        Ls = pi .* exp.(xs)
        ys = [predict(fit, L) for L in Ls]
    end
    plot!(plt, xs, ys; label = label, color = color, linestyle = linestyle, linewidth = 2.2)
    return plt
end

function make_input_plot(observations)
    plt = base_plot(title = "Canonical U/t=2 half-chain input: finite-beta DQMC and ground-state references")
    plot_points!(plt, observations; show_labels = false)
    return plt
end

function make_total_plot(observations, fits)
    plt = base_plot(title = "Total fits: c_eff is a finite-size diagnostic only")
    plot_points!(plt, observations)
    dqmc = principal_fit(fits, "dqmc", "total", ALL_BRANCH)
    without = principal_fit(fits, "reference_without_L26", "total", ALL_BRANCH)
    with = principal_fit(fits, "reference_with_L26_candidate", "total", ALL_BRANCH)
    add_fit_line!(plt, dqmc; label = @sprintf("DQMC WLS c_eff diagnostic=%.3f", dqmc.c_eff_diagnostic), color = COLOR_DQMC)
    add_fit_line!(plt, without; label = @sprintf("Reference without L26 OLS c_eff diagnostic=%.3f", without.c_eff_diagnostic), color = COLOR_WITHOUT, linestyle = :dash)
    add_fit_line!(plt, with; label = @sprintf("Reference with provisional L26 OLS c_eff diagnostic=%.3f", with.c_eff_diagnostic), color = COLOR_WITH, linestyle = :dot)
    return plt
end

function make_branch_plot(observations, fits)
    panels = Any[]
    for branch in (BRANCH_0, BRANCH_2)
        panel = base_plot(title = "$branch diagnostic fits", size = (700, 720), legend = :topleft)
        plot_points!(panel, observations; branch = branch)
        dqmc = principal_fit(fits, "dqmc", "branch", branch)
        without = principal_fit(fits, "reference_without_L26", "branch", branch)
        with = principal_fit(fits, "reference_with_L26_candidate", "branch", branch)
        add_fit_line!(panel, dqmc; label = @sprintf("DQMC WLS c_eff diagnostic=%.3f", dqmc.c_eff_diagnostic), color = COLOR_DQMC)
        add_fit_line!(panel, without; label = @sprintf("Reference without L26 c_eff diagnostic=%.3f", without.c_eff_diagnostic), color = COLOR_WITHOUT, linestyle = :dash)
        add_fit_line!(panel, with; label = @sprintf("Reference with provisional L26 c_eff diagnostic=%.3f", with.c_eff_diagnostic), color = COLOR_WITH, linestyle = :dot)
        push!(panels, panel)
    end
    return plot(panels...; layout = (1, 2), size = (1500, 760), dpi = 180,
        plot_title = "L mod 4 branch fits; c_eff is diagnostic, not an asymptotic determination")
end

function make_oscillation_plot(observations, fits)
    plt = base_plot(title = "Oscillatory fits (fixed p=1 and p=2): finite-size diagnostics only", legend = :outerright)
    plot_points!(plt, observations)
    configurations = [
        ("dqmc", COLOR_DQMC, "DQMC"),
        ("reference_without_L26", COLOR_WITHOUT, "Reference without L26"),
        ("reference_with_L26_candidate", COLOR_WITH, "Reference with provisional L26"),
    ]
    for (dataset, color, label_prefix) in configurations, p in (1, 2)
        fit = principal_fit(fits, dataset, "oscillation", ALL_BRANCH; p = p)
        add_fit_line!(plt, fit; label = @sprintf("%s p=%d c_eff diagnostic=%.3f", label_prefix, p, fit.c_eff_diagnostic),
            color = color, linestyle = p == 1 ? :solid : :dash, discrete = true)
    end
    return plt
end

const FIGURE_STEMS = [
    "oneD_hubbard_fss_fit_input",
    "oneD_hubbard_fss_total_fit",
    "oneD_hubbard_fss_branch_fit",
    "oneD_hubbard_fss_oscillation_fit",
]

function canonicalize_svg_clip_ids(path::String)
    svg = read(path, String)
    source_ids = unique([String(match.captures[1]) for match in eachmatch(r"id=\"(clip[0-9]+)\"", svg)])
    mapping = Dict(source_id => "clip" * lpad(string(index - 1), 3, '0') for (index, source_id) in enumerate(source_ids))
    canonical = replace(svg, r"clip[0-9]+" => token -> mapping[String(token)])
    open(path, "w") do io
        write(io, canonical)
    end
    return nothing
end

function assert_outputs_absent(options::Options)
    ispath(options.summary_output) && error("Refusing to overwrite existing summary output: $(options.summary_output)")
    for stem in FIGURE_STEMS, extension in ("png", "svg")
        target = joinpath(options.figure_dir, "$stem.$extension")
        ispath(target) && error("Refusing to overwrite existing figure output: $target")
    end
end

function render_figures(options::Options, observations, fits)
    plots = [
        make_input_plot(observations),
        make_total_plot(observations, fits),
        make_branch_plot(observations, fits),
        make_oscillation_plot(observations, fits),
    ]
    for (stem, plot_object) in zip(FIGURE_STEMS, plots)
        svg_path = joinpath(options.figure_dir, "$stem.svg")
        savefig(plot_object, svg_path)
        canonicalize_svg_clip_ids(svg_path)
        savefig(plot_object, joinpath(options.figure_dir, "$stem.png"))
    end
end

function validate_principal_count_gates(fits::Vector{FitResult})
    dqmc_points = principal_fit(fits, "dqmc", "total", ALL_BRANCH).n_points
    dqmc_expected = dqmc_points == 8 ?
        (8, Dict(BRANCH_0 => 4, BRANCH_2 => 4), 6, 5, Dict(BRANCH_0 => 2, BRANCH_2 => 2)) :
        dqmc_points == 10 ?
        (10, Dict(BRANCH_0 => 5, BRANCH_2 => 5), 8, 7, Dict(BRANCH_0 => 3, BRANCH_2 => 3)) :
        error("Unexpected principal DQMC point count: $dqmc_points")
    expected = Dict(
        "dqmc" => dqmc_expected,
        "reference_without_L26" => (11, Dict(BRANCH_0 => 6, BRANCH_2 => 5), 9, 8, Dict(BRANCH_0 => 4, BRANCH_2 => 3)),
        "reference_with_L26_candidate" => (12, Dict(BRANCH_0 => 6, BRANCH_2 => 6), 10, 9, Dict(BRANCH_0 => 4, BRANCH_2 => 4)),
    )
    for dataset in DATASET_ORDER
        total_points, branch_counts, total_dof, osc_dof, branch_dof = expected[dataset]
        total = principal_fit(fits, dataset, "total", ALL_BRANCH)
        total.n_points == total_points && total.dof == total_dof || error("Principal total count gate failed for $dataset")
        for p in (1, 2)
            osc = principal_fit(fits, dataset, "oscillation", ALL_BRANCH; p = p)
            osc.n_points == total_points && osc.dof == osc_dof || error("Principal oscillation count gate failed for $dataset p=$p")
        end
        for branch in (BRANCH_0, BRANCH_2)
            fit = principal_fit(fits, dataset, "branch", branch)
            fit.n_points == branch_counts[branch] && fit.dof == branch_dof[branch] || error("Principal branch count gate failed for $dataset/$branch")
        end
    end
end

function print_principal_summary(fits::Vector{FitResult}, input_hash::String, options::Options)
    observation_count = principal_fit(fits, "dqmc", "total", ALL_BRANCH).n_points +
        sum(length, values(EXPECTED_REFERENCE_L))
    println("Validated $observation_count canonical observations; input SHA-256: $input_hash")
    println("Generated $(length(fits)) fit records (principal + Lmin + windows + leave-one-out)")
    println("DQMC covariance: independent absolute 1sigma, Sigma=diag(sigma_i^2), no reduced-chi2 rescaling")
    println("Reference errors: nominal_ols_se from regression scatter; not physical/precision uncertainties")
    println("c_eff is a finite-size diagnostic only, not an asymptotic central-charge determination")
    println()
    println("Principal total diagnostics:")
    for dataset in DATASET_ORDER
        fit = principal_fit(fits, dataset, "total", ALL_BRANCH)
        @printf("  %-32s c_eff=% .9f  %s=% .9g  cond_norm=%.3g\n",
            dataset, fit.c_eff_diagnostic, fit.statistic_name, fit.statistic_value,
            fit.condition_number_normalized)
    end
    println("Summary: ", options.summary_output)
    println("Figures: ", options.figure_dir)
end

function main(args::Vector{String})
    options = parse_options(args)
    options === nothing && return 0
    assert_outputs_absent(options)
    raw_rows = read_csv_dicts(options.input)
    validate_input_rows(raw_rows)
    observations = to_observations(raw_rows)
    input_hash = sha256_file(options.input)
    fits = run_all_fits(observations, input_hash)
    validate_principal_count_gates(fits)

    mkpath(dirname(options.summary_output))
    mkpath(options.figure_dir)
    render_figures(options, observations, fits)
    write_summary(options.summary_output, fits)
    print_principal_summary(fits, input_hash, options)
    return 0
end

try
    exit(main(ARGS))
catch error_object
    println(stderr, "ERROR: ", sprint(showerror, error_object))
    exit(2)
end
