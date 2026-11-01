using LinearAlgebra
using Printf

# L, beta, nA, S2_DQMC, S2_err, S2_ED
rows = [
    (4,  16, 2, 1.113651203984, 0.009975514074400, 0.8803580896857),
    (6,  24, 3, 1.533082262570, 0.002249188431683, 1.535020904893),
    (8,  32, 4, 1.318734225375, 0.007790202477854, 1.309407045626),
    (10, 40, 5, 1.716441955736, 0.003040377094000, 1.716551661433),
    (12, 48, 6, 1.561005400404, 0.007601610411883, 1.542354369358),
]

function weighted_fit(data, label)
    Lvals = Float64[r[1] for r in data]
    y = Float64[r[4] for r in data]
    sigma = Float64[r[5] for r in data]
    x = log.(Lvals)

    X = hcat(x, ones(length(x)))
    W = Diagonal(1.0 ./ sigma.^2)

    cov = inv(transpose(X) * W * X)
    coeff = cov * (transpose(X) * W * y)

    A = coeff[1]
    B = coeff[2]
    Aerr = sqrt(cov[1,1])
    Berr = sqrt(cov[2,2])

    yfit = X * coeff
    chi2 = sum(((y .- yfit) ./ sigma).^2)
    dof = length(y) - 2

    println()
    println(label)
    println("=" ^ length(label))
    @printf("A = %.12f +/- %.12f\n", A, Aerr)
    @printf("B = %.12f +/- %.12f\n", B, Berr)
    @printf("chi2 = %.6f\n", chi2)
    println("dof  = $dof")
    if dof > 0
        @printf("chi2/dof = %.6f\n", chi2 / dof)
    end

    println()
    println("Rows:")
    for (r, yy, yf) in zip(data, y, yfit)
        @printf(
            "L=%2d  logL=%.6f  S2=%.12f  fit=%.12f  residual=%+.12f\n",
            r[1], log(r[1]), yy, yf, yy - yf
        )
    end
end

weighted_fit(rows, "Fit A: include L4, L6, L8, L10, L12")
weighted_fit(filter(r -> r[1] != 4, rows), "Fit B: exclude L4 diagnostic point")
