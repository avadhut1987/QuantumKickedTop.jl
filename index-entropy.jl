using Distributed
addprocs()

@everywhere begin
    using QuantumKickedTop
    using QuantumKickedTop.PhiStates
    using QuantumKickedTop.QuantumUtils
    using QuantumKickedTop.HeatmapEntropy: vn_entropy_grid, density_plot
    using QuantumKickedTop.FloquetSystem: floquet
end

# ========== PARAMETERS ==========
const res = 200
const itr = 1000
θs = LinRange(0.0, π, res)
ϕs = LinRange(-π, π, res)
space = (θs=θs, ϕs=ϕs)

# ========== MAIN LOOP ==========
s = 50.5
Q = 1

# krlist = [1.0]
krlist = [0.99985] .* (s * π / 2)

p = π/2
for kr in krlist
    kt = 1.0

    k = kr + kt
    kp = kr - kt

    params = (p=p, k=k, kp=kp)
    println("▶ Running case: kr = $kr, kt = $kt")


    U = floquet(s, p, k, kp)
    entropy_data, diag = vn_entropy_grid(s, U, Q, itr, space)
    println("Diagnostics: $diag")

    density_plot(entropy_data, Q, space, s, params)
end

# p = 2π
# k = 0.5
# kplist = [3.5, 4.0, 4.5]

# for kp in kplist
#     params = (p=p, k=k, kp=kp)
#     println("▶ Running case: k = $k, kp = $kp")

#     U = floquet(s, p, k, kp)
#     entropy_data, diag = vn_entropy_grid(s, U, Q, itr, space)
#     println("Diagnostics: $diag")

#     density_plot(entropy_data, Q, space, s, params)
# end
