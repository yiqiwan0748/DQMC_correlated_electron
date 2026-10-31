using LinearAlgebra
using Printf

struct Point
    L::Int
    nA::Int
    S2::Float64
    err::Float64
end

# U=2, beta=4L, half-chain, Ninc=16
# L16 has no ED benchmark, but passed internal diagnostics.
points = Point[
    Point(4,  2, 1.113651203984, 0.009975514074400),
    Point(6,  3, 1.533082262570, 0.002249188431683),
    Point(8,  4, 1.318734225375, 0.007790202477854),
    Point(10, 5, 1.716441955736, 0.003040377094000),
    Point(12, 6, 1.561005400404, 0.007601610411883),
    Point(14, 7, 1.843169617829, 0.003662701883375),
    Point(16, 8, 1.705923245750, 0.007900417155705),
]

function weighted_fit(branch)
    x = log.([p.L for p in branch])
    y = [p.S2 for p in branch]
    σ = [p.err for p in branch]

    X = hcat(x, ones(length(x)))
    W = Diagonal(1.0 ./ σ.^2)

    cov = inv(transpose(X) * W * X)
    β = cov * transpose(X) * W * y

    yfit = X * β
    χ2 = sum(((y .- yfit) ./ σ).^2)
    dof = length(y) - 2

    A = β[1]
    B = β[2]
    Aerr = sqrt(cov[1, 1])
    Berr = sqrt(cov[2, 2])

    return A, B, Aerr, Berr, χ2, dof, yfit
end

function write_branch(io, name, branch)
    println(io, name)
    println(io, "L    nA    S2_DQMC              err")
    for p in branch
        @printf(io, "%-4d %-5d %.12f    %.12f\n", p.L, p.nA, p.S2, p.err)
    end

    A, B, Aerr, Berr, χ2, dof, yfit = weighted_fit(branch)

    println(io)
    @printf(io, "Fit form: S2(L) = A log(L) + B\n")
    @printf(io, "A = %.12f +/- %.12f\n", A, Aerr)
    @printf(io, "B = %.12f +/- %.12f\n", B, Berr)
    @printf(io, "chi2 = %.12f\n", χ2)
    @printf(io, "dof = %d\n", dof)
    if dof > 0
        @printf(io, "chi2/dof = %.12f\n", χ2 / dof)
    else
        @printf(io, "chi2/dof = NaN\n")
    end
    @printf(io, "CFT expected A = c/4 = 0.25 for c=1\n")
    @printf(io, "A - 0.25 = %.12f\n", A - 0.25)

    println(io)
    println(io, "Residuals:")
    println(io, "L    residual/err")
    for (p, yf) in zip(branch, yfit)
        @printf(io, "%-4d %.12f\n", p.L, (p.S2 - yf) / p.err)
    end
end

even_branch = [p for p in points if iseven(p.nA)]
odd_branch  = [p for p in points if isodd(p.nA)]

out = "results/processed/FSS_branch_fit_L4_to_L16.txt"
open(out, "w") do io
    println(io, "U=2 half-chain finite-size branch fit")
    println(io, "Data: beta=4L, Ninc=16")
    println(io, "Fit variable: log(L)")
    println(io)

    write_branch(io, "even nA branch: L = 4, 8, 12, 16", even_branch)
    println(io)
    println(io, repeat("=", 72))
    println(io)
    write_branch(io, "odd nA branch: L = 6, 10, 14", odd_branch)
end

println("Wrote ", out)