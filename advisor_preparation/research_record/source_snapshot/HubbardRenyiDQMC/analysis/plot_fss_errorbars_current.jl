using Printf
using LinearAlgebra
using Plots

# Completed FSS data: L, S2_DQMC, S2_error
L = [4, 6, 8, 10, 12, 14]
S2 = [
    1.113651203984,
    1.533082262570,
    1.318734225375,
    1.716441955736,
    1.561005400404,
    1.843169617829,
]
S2err = [
    0.009975514074400,
    0.002249188431683,
    0.007790202477854,
    0.003040377094000,
    0.007601610411883,
    0.003662701883375,
]

x = log.(L)

function weighted_line_fit(x, y, yerr)
    X = hcat(x, ones(length(x)))
    W = Diagonal(1.0 ./ yerr.^2)
    cov = inv(transpose(X) * W * X)
    coeff = cov * (transpose(X) * W * y)
    A = coeff[1]
    B = coeff[2]
    Aerr = sqrt(cov[1, 1])
    Berr = sqrt(cov[2, 2])
    resid = y .- (A .* x .+ B)
    chi2 = sum((resid ./ yerr).^2)
    dof = length(x) - 2
    return A, B, Aerr, Berr, resid, chi2, dof
end

A, B, Aerr, Berr, resid, chi2, dof = weighted_line_fit(x, S2, S2err)

xfit = range(minimum(x) - 0.05, maximum(x) + 0.05, length=200)
yfit = A .* xfit .+ B

figdir = "results/figures"
mkpath(figdir)

# Main S2 plot
p1 = plot(
    xfit,
    yfit,
    label = @sprintf("weighted fit: A=%.4f", A),
    xlabel = "log(L)",
    ylabel = "S2(L/2)",
    title = "FSS data with Monte Carlo error bars, L4 to L14",
    linewidth = 2,
)

idx_lower = findall(l -> mod(l, 4) == 0, L)
idx_upper = findall(l -> mod(l, 4) == 2, L)

scatter!(
    p1,
    x[idx_lower],
    S2[idx_lower],
    yerror = S2err[idx_lower],
    label = "DQMC, L mod 4 = 0",
    markersize = 5,
)

scatter!(
    p1,
    x[idx_upper],
    S2[idx_upper],
    yerror = S2err[idx_upper],
    label = "DQMC, L mod 4 = 2",
    markersize = 5,
)

for i in eachindex(L)
    annotate!(p1, x[i], S2[i], text("L=$(L[i])", 8, :left))
end

savefig(p1, joinpath(figdir, "fss_S2_errorbars_completed_L4_to_L14.png"))

# Residual plot
p2 = scatter(
    x[idx_lower],
    resid[idx_lower],
    yerror = S2err[idx_lower],
    label = "L mod 4 = 0",
    xlabel = "log(L)",
    ylabel = "S2 - weighted fit",
    title = "FSS residuals with Monte Carlo error bars, L4 to L14",
    markersize = 5,
)

scatter!(
    p2,
    x[idx_upper],
    resid[idx_upper],
    yerror = S2err[idx_upper],
    label = "L mod 4 = 2",
    markersize = 5,
)

hline!(p2, [0.0], label = "zero")

for i in eachindex(L)
    annotate!(p2, x[i], resid[i], text("L=$(L[i])", 8, :left))
end

savefig(p2, joinpath(figdir, "fss_residuals_errorbars_completed_L4_to_L14.png"))

# Error-bar-size plot
p3 = bar(
    string.(L),
    S2err,
    label = false,
    xlabel = "L",
    ylabel = "S2 error bar",
    title = "Monte Carlo error size by system size, L4 to L14",
)

savefig(p3, joinpath(figdir, "fss_errorbar_size_completed_L4_to_L14.png"))

println("Saved figures:")
println(joinpath(figdir, "fss_S2_errorbars_completed_L4_to_L14.png"))
println(joinpath(figdir, "fss_residuals_errorbars_completed_L4_to_L14.png"))
println(joinpath(figdir, "fss_errorbar_size_completed_L4_to_L14.png"))

println()
println("Weighted fit to all completed points L4 to L14:")
@printf("A = %.12f +/- %.12f\n", A, Aerr)
@printf("B = %.12f +/- %.12f\n", B, Berr)
@printf("chi2/dof = %.6f / %d = %.6f\n", chi2, dof, chi2 / dof)

println()
println("Data table:")
for i in eachindex(L)
    @printf("L=%2d  logL=%.6f  S2=%.12f +/- %.12f  residual=%+.12f\n",
            L[i], x[i], S2[i], S2err[i], resid[i])
end
