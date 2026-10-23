using Plots
using Printf

figdir = "results/figures"
mkpath(figdir)

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

L = [r[1] for r in data]
nA = [r[2] for r in data]
x = log.(L)
S2 = [r[3] for r in data]
S2err = [r[4] for r in data]

A_even = 0.407393316576
B_even = 0.519740650716

A_odd = 0.364246342389
B_odd = 0.879957216099

idx_lower = findall(l -> mod(l, 4) == 0, L)  # even nA
idx_upper = findall(l -> mod(l, 4) == 2, L)  # odd nA

xfit = range(log(4), log(14), length=200)

# Dummy/reference first series keeps the same default color-cycle convention
# as plot_fss_errorbars_current.jl:
#   series 1 = reference line
#   series 2 = lower branch, L mod 4 = 0
#   series 3 = upper branch, L mod 4 = 2
p1 = plot(
    xfit,
    0.25 .* xfit .+ 0.95,
    label = "reference slope 0.25",
    xlabel = "log(L)",
    ylabel = "S2(L/2)",
    title = "FSS branch fits, L4 to L14",
    linewidth = 2,
    linestyle = :dash,
)

plot!(
    p1,
    xfit,
    A_even .* xfit .+ B_even,
    label = @sprintf("L mod 4 = 0 fit: A=%.3f", A_even),
    linewidth = 2,
)

plot!(
    p1,
    xfit,
    A_odd .* xfit .+ B_odd,
    label = @sprintf("L mod 4 = 2 fit: A=%.3f", A_odd),
    linewidth = 2,
)

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

savefig(p1, joinpath(figdir, "fss_branch_fits_L4_to_L14.png"))

println("Saved figure:")
println(joinpath(figdir, "fss_branch_fits_L4_to_L14.png"))
