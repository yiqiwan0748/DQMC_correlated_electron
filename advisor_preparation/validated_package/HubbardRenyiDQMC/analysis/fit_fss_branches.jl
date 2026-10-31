using LinearAlgebra
using Printf

# Completed U=2, beta=4L, nA=L/2, Ninc=16 FSS data
# columns: L, nA, S2_DQMC, error
data = [
    (4,  2, 1.113651203984, 0.009975514074400),
    (6,  3, 1.533082262570, 0.002249188431683),
    (8,  4, 1.318734225375, 0.007790202477854),
    (10, 5, 1.716441955736, 0.003040377094000),
    (12, 6, 1.561005400404, 0.007601610411883),
    (14, 7, 1.843169617829, 0.003662701883375),
]

function weighted_log_fit(rows)
    L  = Float64[r[1] for r in rows]
    y  = Float64[r[3] for r in rows]
    dy = Float64[r[4] for r in rows]

    x = log.(L)
    X = hcat(x, ones(length(x)))
    W = Diagonal(1.0 ./ dy.^2)

    cov = inv(transpose(X) * W * X)
    p = cov * (transpose(X) * W * y)

    A, B = p
    dA, dB = sqrt.(diag(cov))

    residuals = y - X * p
    chi2 = sum((residuals ./ dy).^2)
    dof = length(y) - 2

    return A, B, dA, dB, chi2, dof, residuals
end

even_nA = [r for r in data if iseven(r[2])]
odd_nA  = [r for r in data if isodd(r[2])]

println()
println("FSS branch fits for U=2, beta=4L, nA=L/2, Ninc=16")
println("Fit form: S2(L/2) = A log(L) + B")
println("Theory for S2 with c=1: A = c/4 = 0.25")
println()

for (label, rows) in [("even nA", even_nA), ("odd nA", odd_nA)]
    A, B, dA, dB, chi2, dof, residuals = weighted_log_fit(rows)

    println(label)
    println("points:")
    for r in rows
        @printf("  L=%2d  nA=%2d  S2=%.12f +/- %.12f\n", r[1], r[2], r[3], r[4])
    end

    @printf("A = %.12f +/- %.12f\n", A, dA)
    @printf("B = %.12f +/- %.12f\n", B, dB)
    @printf("chi2/dof = %.6f / %d = %.6f\n", chi2, dof, chi2/dof)
    @printf("A - 0.25 = %.12f\n", A - 0.25)
    println()
end
