#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_SUMMARY = joinpath(
    PROJECT_ROOT,
    "results",
    "processed",
    "BLOCKSIZE_U2_beta40_Ninc16_summary.csv",
)
const EXPECTED_HEADER = [
    "nA",
    "S2_DQMC",
    "S2_error",
    "S2_ED",
    "DeltaS2",
    "z",
    "status",
]

function usage(io::IO = stdout)
    println(io, "Usage:")
    println(io, "  julia analysis/plot_L10_U2_beta40_blocksize.jl [--input-summary PATH] --output PATH")
    println(io)
    println(io, "Options:")
    println(io, "  --input-summary PATH  Tracked L10 block-size summary CSV")
    println(io, "                        (default: results/processed/BLOCKSIZE_U2_beta40_Ninc16_summary.csv)")
    println(io, "  --output PATH         New SVG output path; existing paths are refused")
    println(io, "  --help                Show this help and exit")
end

function parse_cli(args::Vector{String})
    if any(==("--help"), args)
        args == ["--help"] || error("--help cannot be combined with other arguments")
        return (help = true, input = "", output = "")
    end

    input = DEFAULT_SUMMARY
    output = nothing
    seen_input = false
    seen_output = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--input-summary"
            seen_input && error("--input-summary may be provided only once")
            i == length(args) && error("--input-summary requires a path")
            input = args[i + 1]
            isempty(input) && error("--input-summary requires a nonempty path")
            seen_input = true
            i += 2
        elseif arg == "--output"
            seen_output && error("--output may be provided only once")
            i == length(args) && error("--output requires a path")
            output = args[i + 1]
            isempty(output) && error("--output requires a nonempty path")
            seen_output = true
            i += 2
        else
            error("Unknown argument: $arg")
        end
    end

    output === nothing && error("--output is required")
    input_path = abspath(input)
    output_path = abspath(output)
    isfile(input_path) || error("Input summary is not a file")
    (ispath(output_path) || islink(output_path)) && error("Output path already exists")
    isdir(dirname(output_path)) || error("Output directory does not exist")
    lowercase(splitext(output_path)[2]) == ".svg" || error("Output path must end in .svg")
    return (help = false, input = input_path, output = output_path)
end

function parse_finite_float(token::AbstractString, field::String, line_number::Int)
    value = tryparse(Float64, token)
    value === nothing && error("Invalid $field at CSV line $line_number: $token")
    isfinite(value) || error("Non-finite $field at CSV line $line_number: $token")
    return value
end

function read_summary(path::String)
    lines = readlines(path)
    isempty(lines) && error("Input summary is empty")

    header = strip.(split(lines[1], ","; keepempty = true))
    header == EXPECTED_HEADER || error(
        "Unexpected CSV header; expected " * join(EXPECTED_HEADER, ","),
    )

    rows = NamedTuple[]
    for (offset, raw_line) in enumerate(lines[2:end])
        line_number = offset + 1
        isempty(strip(raw_line)) && error("Blank CSV line at line $line_number")
        occursin('"', raw_line) && error("Quoted CSV fields are unsupported at line $line_number")
        fields = strip.(split(raw_line, ","; keepempty = true))
        length(fields) == length(EXPECTED_HEADER) || error(
            "Malformed CSV line $line_number: expected $(length(EXPECTED_HEADER)) fields",
        )

        nA = tryparse(Int, fields[1])
        nA === nothing && error("Invalid nA at CSV line $line_number: $(fields[1])")
        s2_dqmc = parse_finite_float(fields[2], "S2_DQMC", line_number)
        s2_error = parse_finite_float(fields[3], "S2_error", line_number)
        s2_ed = parse_finite_float(fields[4], "S2_ED", line_number)
        delta_stored = parse_finite_float(fields[5], "DeltaS2", line_number)
        z_stored = parse_finite_float(fields[6], "z", line_number)
        isempty(fields[7]) && error("Empty status at CSV line $line_number")
        s2_error > 0 || error("S2_error must be positive at CSV line $line_number")

        delta_expected = s2_dqmc - s2_ed
        isapprox(delta_stored, delta_expected; atol = 1.1e-12, rtol = 1e-12) || error(
            "DeltaS2 is inconsistent at CSV line $line_number",
        )
        z_expected = delta_expected / s2_error
        isapprox(z_stored, z_expected; atol = 5e-4, rtol = 1e-12) || error(
            "z is inconsistent at CSV line $line_number",
        )

        push!(rows, (
            nA = nA,
            S2_DQMC = s2_dqmc,
            S2_error = s2_error,
            S2_ED = s2_ed,
        ))
    end

    sort!(rows; by = row -> row.nA)
    [row.nA for row in rows] == collect(1:5) || error(
        "Input must contain exactly one row for each nA=1:5",
    )
    return rows
end

svg_escape(value) = replace(
    string(value),
    "&" => "&amp;",
    "<" => "&lt;",
    ">" => "&gt;",
    "\"" => "&quot;",
)

function svg_text(x, y, text; anchor = "middle", size = 14, weight = "normal", extra = "")
    return @sprintf(
        "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"%s\" font-family=\"Arial, sans-serif\" font-size=\"%d\" font-weight=\"%s\" %s>%s</text>",
        x,
        y,
        anchor,
        size,
        weight,
        extra,
        svg_escape(text),
    )
end

diamond_points(x, y, radius) = @sprintf(
    "%.2f,%.2f %.2f,%.2f %.2f,%.2f %.2f,%.2f",
    x,
    y - radius,
    x + radius,
    y,
    x,
    y + radius,
    x - radius,
    y,
)

function make_svg(rows)
    width = 920.0
    height = 620.0
    left = 92.0
    right = 32.0
    top = 72.0
    bottom = 88.0
    plot_width = width - left - right
    plot_height = height - top - bottom

    y_values = vcat(
        [row.S2_DQMC - row.S2_error for row in rows],
        [row.S2_DQMC + row.S2_error for row in rows],
        [row.S2_ED for row in rows],
    )
    ymin_tick = floor(Int, (minimum(y_values) - 0.04) * 10)
    ymax_tick = ceil(Int, (maximum(y_values) + 0.04) * 10)
    ymin = ymin_tick / 10
    ymax = ymax_tick / 10
    xmin = 1.0 - 0.35
    xmax = 5.0 + 0.35

    sx(x) = left + (x - xmin) / (xmax - xmin) * plot_width
    sy(y) = top + (ymax - y) / (ymax - ymin) * plot_height

    out = IOBuffer()
    println(out, @sprintf(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" viewBox=\"0 0 %.0f %.0f\" role=\"img\" aria-labelledby=\"plot-title plot-description\">",
        width,
        height,
        width,
        height,
    ))
    println(out, "<title id=\"plot-title\">L=10 finite-beta DQMC versus ground-state ED</title>")
    println(out, "<desc id=\"plot-description\">Renyi-2 entropy for subsystem sizes nA=1 through 5; DQMC error bars are one Monte Carlo standard error.</desc>")
    println(out, "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>")
    println(out, svg_text(width / 2, 34, "L=10, U/t=2, beta=40 block-size Renyi-2 benchmark"; size = 22, weight = "700"))
    println(out, svg_text(width / 2, 57, "finite-beta DQMC (Ninc=16) versus ground-state ED"; size = 14))
    println(out, @sprintf(
        "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"#ffffff\" stroke=\"#333333\" stroke-width=\"1\"/>",
        left,
        top,
        plot_width,
        plot_height,
    ))

    for tick in ymin_tick:ymax_tick
        y = tick / 10
        py = sy(y)
        println(out, @sprintf("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>", left, py, left + plot_width, py))
        println(out, @sprintf("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#333333\" stroke-width=\"1\"/>", left - 5, py, left, py))
        println(out, svg_text(left - 10, py + 5, @sprintf("%.1f", y); anchor = "end", size = 12))
    end

    for x in 1:5
        px = sx(x)
        println(out, @sprintf("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>", px, top, px, top + plot_height))
        println(out, @sprintf("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#333333\" stroke-width=\"1\"/>", px, top + plot_height, px, top + plot_height + 5))
        println(out, svg_text(px, top + plot_height + 24, string(x); size = 13))
    end

    println(out, svg_text(left + plot_width / 2, height - 28, "subsystem size nA"; size = 15))
    println(out, svg_text(23, top + plot_height / 2, "S2"; size = 15, extra = @sprintf("transform=\"rotate(-90 %.2f %.2f)\"", 23.0, top + plot_height / 2)))

    dqmc_polyline = join([@sprintf("%.2f,%.2f", sx(row.nA), sy(row.S2_DQMC)) for row in rows], " ")
    println(out, "<polyline class=\"dqmc-series\" points=\"$dqmc_polyline\" fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"2.2\"/>")
    for row in rows
        x = sx(row.nA)
        y = sy(row.S2_DQMC)
        ylow = sy(row.S2_DQMC - row.S2_error)
        yhigh = sy(row.S2_DQMC + row.S2_error)
        println(out, @sprintf("<g class=\"dqmc-errorbar\" data-nA=\"%d\" data-sigma=\"%.12f\">", row.nA, row.S2_error))
        println(out, @sprintf("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#1f77b4\" stroke-width=\"1.8\"/>", x, yhigh, x, ylow))
        println(out, @sprintf("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#1f77b4\" stroke-width=\"1.8\"/>", x - 8, yhigh, x + 8, yhigh))
        println(out, @sprintf("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#1f77b4\" stroke-width=\"1.8\"/>", x - 8, ylow, x + 8, ylow))
        println(out, "</g>")
        println(out, @sprintf("<circle class=\"dqmc-point\" data-nA=\"%d\" data-S2=\"%.12f\" cx=\"%.2f\" cy=\"%.2f\" r=\"5.5\" fill=\"#1f77b4\" stroke=\"#ffffff\" stroke-width=\"1.5\"/>", row.nA, row.S2_DQMC, x, y))
    end

    for row in rows
        x = sx(row.nA)
        y = sy(row.S2_ED)
        println(out, @sprintf("<polygon class=\"ed-point\" data-nA=\"%d\" data-S2=\"%.12f\" points=\"%s\" fill=\"#d62728\" stroke=\"#ffffff\" stroke-width=\"1.2\"/>", row.nA, row.S2_ED, diamond_points(x, y, 7.0)))
    end

    legend_x = left + 20
    legend_y = top + 20
    println(out, @sprintf("<rect x=\"%.2f\" y=\"%.2f\" width=\"264\" height=\"74\" fill=\"#ffffff\" stroke=\"#cccccc\"/>", legend_x, legend_y))
    println(out, @sprintf("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#1f77b4\" stroke-width=\"2\"/>", legend_x + 14, legend_y + 22, legend_x + 46, legend_y + 22))
    println(out, @sprintf("<circle cx=\"%.2f\" cy=\"%.2f\" r=\"5\" fill=\"#1f77b4\"/>", legend_x + 30, legend_y + 22))
    println(out, svg_text(legend_x + 56, legend_y + 27, "DQMC beta=40 +/- 1 MC sigma"; anchor = "start", size = 13))
    println(out, @sprintf("<polygon points=\"%s\" fill=\"#d62728\"/>", diamond_points(legend_x + 30, legend_y + 49, 7.0)))
    println(out, svg_text(legend_x + 56, legend_y + 54, "ground-state ED"; anchor = "start", size = 13))
    println(out, "</svg>")
    return String(take!(out))
end

function write_exclusive(path::String, content::String)
    flags = Base.Filesystem.JL_O_WRONLY |
            Base.Filesystem.JL_O_CREAT |
            Base.Filesystem.JL_O_EXCL |
            Base.Filesystem.JL_O_CLOEXEC
    file = try
        Base.Filesystem.open(path, flags, 0o644)
    catch err
        (ispath(path) || islink(path)) && error("Output path already exists")
        rethrow(err)
    end
    try
        write(file, content)
    finally
        close(file)
    end
    chmod(path, 0o644)
end

function main(args::Vector{String})
    options = parse_cli(args)
    if options.help
        usage()
        return 0
    end
    rows = read_summary(options.input)
    write_exclusive(options.output, make_svg(rows))
    println("Wrote deterministic SVG.")
    return 0
end

try
    exit(main(ARGS))
catch err
    println(stderr, "ERROR: ", sprint(showerror, err))
    exit(1)
end
