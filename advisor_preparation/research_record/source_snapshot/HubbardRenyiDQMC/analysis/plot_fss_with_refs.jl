using Printf
using Plots

figdir = "results/figures"
mkpath(figdir)

# ------------------------------------------------------------
# DQMC beta = 4L data.
# These are finite-beta DQMC results with statistical error bars.
# ------------------------------------------------------------

L_dqmc = [4, 6, 8, 10, 12, 14, 16, 18]

S2_dqmc = [
    1.113651203984,
    1.533082262570,
    1.318734225375,
    1.716441955736,
    1.561005400404,
    1.843169617829,
    1.705923245750,
    1.949714833089,
]

err_dqmc = [
    0.009975514074400,
    0.002249188431683,
    0.007790202477854,
    0.003040377094000,
    0.007601610411883,
    0.003662701883375,
    0.007900417155705,
    0.004306010013,
]

# ------------------------------------------------------------
# ED ground-state references.
# ------------------------------------------------------------

L_ed = [4, 6, 8, 10, 12]

S2_ed = [
    0.8803580896857,
    1.535020904892800,
    1.309407045625562,
    1.716551661433,
    1.542354369357964,
]

# ------------------------------------------------------------
# DMRG/MPS ground-state references.
# L18/L20/L22/L24 are first-pass references because maxlinkdim reached 4800.
# ------------------------------------------------------------

L_dmrg = [14, 16, 18, 20, 22, 24]

S2_dmrg = [
    1.8460982712452576,
    1.6962908755077759,
    1.9442405246094483,
    1.8087666550463910,
    2.0219147921796217,
    1.8960694735263326,
]

# ------------------------------------------------------------
# Work in x = log(L).
# ------------------------------------------------------------

x_dqmc = log.(L_dqmc)
x_ed = log.(L_ed)
x_dmrg = log.(L_dmrg)

idx_dqmc_0 = findall(L -> mod(L, 4) == 0, L_dqmc)
idx_dqmc_2 = findall(L -> mod(L, 4) == 2, L_dqmc)

idx_ed_0 = findall(L -> mod(L, 4) == 0, L_ed)
idx_ed_2 = findall(L -> mod(L, 4) == 2, L_ed)

idx_dmrg_0 = findall(L -> mod(L, 4) == 0, L_dmrg)
idx_dmrg_2 = findall(L -> mod(L, 4) == 2, L_dmrg)

# CFT slope reference A = 1/4, anchored at the L24 DMRG point.
xfit = range(log(4), log(24), length = 200)
A_cft = 0.25
x0 = log(24)
y0 = S2_dmrg[end]
y_cft = y0 .+ A_cft .* (xfit .- x0)

p = plot(
    xfit,
    y_cft,
    label = "CFT slope 1/4 reference",
    linestyle = :dash,
    linewidth = 2,
    xlabel = "L, plotted on log scale",
    ylabel = "S2",
    title = "1D Hubbard U=2 half-chain S2 FSS",
    legend = :topleft,
    size = (950, 640),
    dpi = 300,
)

scatter!(
    p,
    x_dqmc[idx_dqmc_0],
    S2_dqmc[idx_dqmc_0],
    yerror = err_dqmc[idx_dqmc_0],
    label = "DQMC beta=4L, L mod 4 = 0",
    marker = :circle,
    markersize = 5,
)

scatter!(
    p,
    x_dqmc[idx_dqmc_2],
    S2_dqmc[idx_dqmc_2],
    yerror = err_dqmc[idx_dqmc_2],
    label = "DQMC beta=4L, L mod 4 = 2",
    marker = :diamond,
    markersize = 5,
)

scatter!(
    p,
    x_ed[idx_ed_0],
    S2_ed[idx_ed_0],
    label = "ED, L mod 4 = 0",
    marker = :xcross,
    markersize = 7,
)

scatter!(
    p,
    x_ed[idx_ed_2],
    S2_ed[idx_ed_2],
    label = "ED, L mod 4 = 2",
    marker = :cross,
    markersize = 7,
)

scatter!(
    p,
    x_dmrg[idx_dmrg_0],
    S2_dmrg[idx_dmrg_0],
    label = "DMRG/MPS, L mod 4 = 0",
    marker = :square,
    markersize = 6,
)

scatter!(
    p,
    x_dmrg[idx_dmrg_2],
    S2_dmrg[idx_dmrg_2],
    label = "DMRG/MPS, L mod 4 = 2",
    marker = :utriangle,
    markersize = 6,
)

all_L_ticks = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
xticks!(p, log.(all_L_ticks), string.(all_L_ticks))

out = joinpath(figdir, "fss_with_ED_DMRG_refs_L4_to_L24.png")
savefig(p, out)

println("Saved figure:")
println(out)

println()
println("Reference data included:")
println("  DQMC beta=4L: L = ", L_dqmc)
println("  ED:           L = ", L_ed)
println("  DMRG/MPS:     L = ", L_dmrg)
println()
println("Note: L18/L20/L22/L24 DMRG are first-pass references because maxlinkdim reached 4800.")