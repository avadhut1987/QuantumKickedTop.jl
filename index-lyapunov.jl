using Distributed
addprocs(min(7, Sys.CPU_THREADS ÷ 2))  # Use 50-75% of cores


@everywhere begin
    using QuantumKickedTop
    using QuantumKickedTop.ClassicalUtils: get_LLE
    using QuantumKickedTop.HeatmapLLE: lle_grid, density_plot_lle
end

# ========== PARAMETERS ==========
res = 500
time = 1000
θs = collect(range(0.0, π, length=res))
ϕs = collect(range(-π, π, length=res))
space = (θs=θs, ϕs=ϕs)

p = 0.0
kprimeList = [3.9, 3.8, 3.7, 3.6]
k = 0.5

# ========== MAIN LOOP ==========
for kprime in kprimeList
    println("▶ Running case: p = $(p), k = $k, k' = $kprime")
    params = (p=p, k=k, kprime=kprime)

    result = lle_grid(space, params, time, res; max_bad_per_point=3)
    density_plot_lle(result, space, params)
    println("LLE heatmap saved!")
end
