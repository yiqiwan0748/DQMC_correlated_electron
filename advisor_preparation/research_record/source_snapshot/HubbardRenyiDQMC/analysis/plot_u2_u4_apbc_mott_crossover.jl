#!/usr/bin/env julia

using Printf

function usage(io::IO = stdout)
    println(io, "Usage:")
    println(io, "  julia analysis/plot_u2_u4_apbc_mott_crossover.jl --output PATH.svg")
end

function parse_output(args)
    args == ["--help"] && return nothing
    length(args) == 2 && args[1] == "--output" || error("Expected --output PATH.svg")
    output = abspath(args[2])
    lowercase(splitext(output)[2]) == ".svg" || error("Output path must end in .svg")
    isdir(dirname(output)) || error("Output directory does not exist")
    (ispath(output) || islink(output)) && error("Output path already exists")
    return output
end

c_eff(s8, s12) = 4 * (s12 - s8) / log(12 / 8)
c_eff_error(e8, e12) = 4 * hypot(e8, e12) / log(12 / 8)

function make_svg()
    U = [2.0, 4.0]
    ed = [
        c_eff(1.438139939952014, 1.639720441049806),
        c_eff(1.145870324242208, 1.258158791329264),
    ]
    dqmc = [
        c_eff(1.446183297421, 1.657932218556),
        c_eff(1.185789049266, 1.349354013793),
    ]
    dqmc_error = [
        c_eff_error(0.008290348556, 0.010169092474),
        c_eff_error(0.020607329183, 0.031463441265),
    ]

    width, height = 820.0, 540.0
    left, right, top, bottom = 92.0, 34.0, 74.0, 82.0
    plot_width = width - left - right
    plot_height = height - top - bottom
    xmin, xmax = 1.5, 4.5
    ymin, ymax = 0.65, 2.55
    sx(x) = left + (x - xmin) / (xmax - xmin) * plot_width
    sy(y) = top + (ymax - y) / (ymax - ymin) * plot_height

    io = IOBuffer()
    println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"820\" height=\"540\" viewBox=\"0 0 820 540\" role=\"img\" aria-labelledby=\"title desc\">")
    println(io, "<title id=\"title\">Closed-shell APBC Mott crossover from U=2 to U=4</title>")
    println(io, "<desc id=\"desc\">Exact ED effective central charge falls from 1.989 at U=2 to 1.108 at U=4. DQMC screening estimates with error bars show the same direction.</desc>")
    println(io, "<rect width=\"820\" height=\"540\" fill=\"white\"/>")
    println(io, @sprintf("<text x=\"%.1f\" y=\"34\" text-anchor=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"20\" font-weight=\"bold\">Closed-shell APBC Mott crossover</text>", width / 2))
    println(io, @sprintf("<text x=\"%.1f\" y=\"56\" text-anchor=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"13\" fill=\"#444\">two-point slope from L=8 and L=12</text>", width / 2))

    for tick in (0.75, 1.0, 1.5, 2.0, 2.5)
        y = sy(tick)
        println(io, @sprintf("<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#e3e3e3\"/>", left, y, width - right, y))
        println(io, @sprintf("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\" dominant-baseline=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"12\">%.2g</text>", left - 10, y, tick))
    end

    for (value, label) in ((1.0, "c = 1"), (2.0, "c = 2"))
        y = sy(value)
        println(io, @sprintf("<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#777\" stroke-width=\"1.4\" stroke-dasharray=\"7 6\"/>", left, y, width - right, y))
        println(io, @sprintf("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#555\">%s</text>", width - right - 7, y - 7, label))
    end

    println(io, @sprintf("<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"none\" stroke=\"#222\" stroke-width=\"1.2\"/>", left, top, plot_width, plot_height))
    for u in U
        x = sx(u)
        println(io, @sprintf("<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#222\"/>", x, height - bottom, x, height - bottom + 6))
        println(io, @sprintf("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"13\">%.0f</text>", x, height - bottom + 24, u))
    end
    println(io, @sprintf("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"15\">U/t</text>", left + plot_width / 2, height - 26))
    println(io, @sprintf("<text x=\"25\" y=\"%.1f\" text-anchor=\"middle\" transform=\"rotate(-90 25 %.1f)\" font-family=\"Arial, sans-serif\" font-size=\"15\">c_eff diagnostic</text>", top + plot_height / 2, top + plot_height / 2))

    println(io, @sprintf("<polyline points=\"%.1f,%.1f %.1f,%.1f\" fill=\"none\" stroke=\"#0072B2\" stroke-width=\"2.5\"/>", sx(U[1]) - 10, sy(ed[1]), sx(U[2]) - 10, sy(ed[2])))
    for i in eachindex(U)
        x, y = sx(U[i]) - 10, sy(ed[i])
        println(io, @sprintf("<circle cx=\"%.1f\" cy=\"%.1f\" r=\"6\" fill=\"#0072B2\" stroke=\"white\" stroke-width=\"1.5\"/>", x, y))
        anchor = i == 1 ? "end" : "start"
        dx = i == 1 ? -10 : 10
        println(io, @sprintf("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"%s\" font-family=\"Arial, sans-serif\" font-size=\"12\" font-weight=\"bold\" fill=\"#0072B2\">%.3f</text>", x + dx, y - 10, anchor, ed[i]))
    end

    for i in eachindex(U)
        x, y = sx(U[i]) + 10, sy(dqmc[i])
        ytop, ybottom = sy(dqmc[i] + dqmc_error[i]), sy(dqmc[i] - dqmc_error[i])
        println(io, @sprintf("<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#D55E00\" stroke-width=\"2\"/>", x, ytop, x, ybottom))
        for cap_y in (ytop, ybottom)
            println(io, @sprintf("<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#D55E00\" stroke-width=\"2\"/>", x - 6, cap_y, x + 6, cap_y))
        end
        points = @sprintf("%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f", x, y - 7, x + 7, y, x, y + 7, x - 7, y)
        println(io, "<polygon points=\"$points\" fill=\"#D55E00\" stroke=\"white\" stroke-width=\"1.5\"/>")
        println(io, @sprintf("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"start\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#D55E00\">%.2f ± %.2f</text>", x + 12, y + 5, dqmc[i], dqmc_error[i]))
    end

    legend_y = top + 20
    println(io, @sprintf("<circle cx=\"%.1f\" cy=\"%.1f\" r=\"5\" fill=\"#0072B2\"/>", left + 166, legend_y))
    println(io, @sprintf("<text x=\"%.1f\" y=\"%.1f\" dominant-baseline=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"12\">Exact ED</text>", left + 178, legend_y))
    println(io, @sprintf("<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"#D55E00\"/>", left + 304, legend_y - 6, left + 310, legend_y, left + 304, legend_y + 6, left + 298, legend_y))
    println(io, @sprintf("<text x=\"%.1f\" y=\"%.1f\" dominant-baseline=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"12\">Ninc=8 DQMC screening</text>", left + 318, legend_y))
    println(io, "</svg>")
    return String(take!(io))
end

function main(args = ARGS)
    output = parse_output(args)
    if output === nothing
        usage()
        return
    end
    write(output, make_svg())
    println("Wrote ", output)
end

main()
