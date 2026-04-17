using Distributed
addprocs()

@everywhere begin
    using QuantumKickedTop
    using QuantumKickedTop.FloquetSystem: floquet
    using QuantumKickedTop.HeatmapDiscord: discord_grid, density_plot_discord
end

const res = 200
const itr = 1000
θs = LinRange(0, π, res)
ϕs = LinRange(-π, π, res)
space = (θs=θs, ϕs=ϕs)

s = 50.5
krlist = [0.999, 0.99, 0.95, 0.9] .* (s * π / 2)

for kr in krlist
    kt = 0.0
    params = (kr=kr, kt=kt)
    println("kr :", kr)

    unitaryOperator = floquet(s, kr, kt)

    # --- NEW: discord on 2-qubit subsystem (Ns=2, split into A=1 qubit, B=1 qubit)
    discord_data = discord_grid(s, unitaryOperator, Ns=2, Q_A=1, itr, space)
    density_plot_discord(discord_data, Ns=2, Q_A=1, space, s, params)
end
