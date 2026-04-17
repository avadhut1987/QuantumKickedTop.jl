using Distributed
addprocs()

@everywhere begin
    using QuantumKickedTop
    using QuantumKickedTop.PhiStates
    using QuantumKickedTop.QuantumUtils
    using QuantumKickedTop.HeatmapLinearEntropy: linear_entropy_grid, density_plot
    using QuantumKickedTop.FloquetSystem: floquet
end



# ========== PARAMETERS ==========
const res = 200
const itr = 1000
θs = LinRange(0, π, res)
ϕs = LinRange(-π, π, res)
space = (θs=θs, ϕs=ϕs)

# ========== MAIN LOOP ==========
s = 50.5
# krlist = [0.9996, 0.9997, 0.99975, 1.0] .* (s * π / 2)
krlist = [1.0, 1.25, 1.5]

Q = 1

for kr in krlist
    kt = kr
    params = (kr=kr, kt=kt)
    println("▶ Running case: kr = $kr, kt = $kt")

    U = floquet(s, kr, kt)
    entropy_data, diag = linear_entropy_grid(s, U, Q, itr, space)
    println("Diagnostics: $diag")

    density_plot(entropy_data, Q, space, s, params)
end
