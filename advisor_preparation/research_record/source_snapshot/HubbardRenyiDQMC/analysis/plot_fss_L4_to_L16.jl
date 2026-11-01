using LinearAlgebra
using Printf
using Plots

L = [4, 6, 8, 10, 12, 14, 16]
nA = [2, 3, 4, 5, 6, 7, 8]

S2 = [
    1.113651203984,
    1.533082262570,
    1.318734225375,
    1.716441955736,
    1.561005400404,
    1.843169617829,
    1.705923245750,
]

S2err = [
    0.009975514074400,
    0.002249188431683,
    0.007790202477854,
    0.003040377094000,
    0.007601610411883,
    0.003662701883375,
    0.007900417155705,
]

function weighted_line_fit(x, y, err)
    X = hcat(x, ones(length(x)))
    W = Diagonal(1.0 ./ err.^2)
    cov = inv(X' * W * X)
    β = cov * X' * W * y
    yfit = X * β
    chi2 = sum(((y .- yfit) ./ err).^2)
    dof = length(y) - 2
    return β[1], β[2], sqrt(cov[1,1]), sqrt(cov[2,2]), chi2, dof
end

x = log.(L)

idx_0 = findall(l -> mod(l, 4) == 0, L)
idx_2 = findall(l -> mod(l, 4) == 2, L)

A0, B0, dA0, dB0, chi20, dof0 = weighted_line_fit(x[idx_0], S2[idx_0], S2err[idx_0])
A2, B2, dA2, dB2, chi22, dof2 = weighted_line_fit(x[idx_2], S2[idx_2], S2err[idx_2])

xfit = range(log(4), log(16), length=200)

figdir = "results/figures"
mkpath(figdir)

p1 = plot(
    xfit,
    0.25 .* xfit .+ 0.95,
    label = "reference slope 0.25",
    xlabel = "log(L)",
    ylabel = "S2(L/2)",
    title = "FSS fits, L4 to L16",
    linewidth = 2,
    linestyle = :dash,
)

scatter!(
    p1,
    x[idx_0],
    S2[idx_0],
    yerror = S2err[idx_0],
    label = "DQMC, L mod 4 = 0",
    markersize = 5,
)

scatter!(
    p1,
    x[idx_2],
    S2[idx_2],
    yerror = S2err[idx_2],
    label = "DQMC, L mod 4 = 2",
    markersize = 5,
)

plot!(
    p1,
    xfit,
    A0 .* xfit .+ B0,
    label = @sprintf("fit, L mod 4 = 0: A = %.3f", A0),
    linewidth = 2,
)

plot!(
    p1,
    xfit,
    A2 .* xfit .+ B2,
    label = @sprintf("fit, L mod 4 = 2: A = %.3f", A2),
    linewidth = 2,
)

savefig(p1, joinpath(figdir, "fss_fits_L4_to_L16.png"))

println("Fit results:")
@printf("L mod 4 = 0: A = %.12f +/- %.12f, chi2/dof = %.12f\n", A0, dA0, chi20 / dof0)
@printf("L mod 4 = 2: A = %.12f +/- %.12f, chi2/dof = %.12f\n", A2, dA2, chi22 / dof2)

println()
println("Saved figure:")
println(joinpath(figdir, "fss_fits_L4_to_L16.png"))