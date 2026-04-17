using QuantumKickedTop
using QuantumKickedTop.QuantumUtils
using QuantumKickedTop.FloquetSystem
using QuantumKickedTop.LinearEntropyVsKR: long_time_avg_entropy_vs_kr
using Plots, NPZ

# PARAMETERS
krlist = range(0.5, 3.0, length=41)    # edit as needed
s_list = [1.5, 10.5, 50.5]
p = π/2                                 # edit as needed
rdm_dim = 1
itr = 1000                              # test mode (speed); set to 1000+ for final

theta0, phi0 = π / 2, -π / 2
results = Dict{Float64,Vector{Float64}}()

for s in s_list
    ent = long_time_avg_entropy_vs_kr(s, p, collect(krlist), rdm_dim, itr;
        theta0=theta0, phi0=phi0)
    results[s] = ent
    npzwrite("npy/entropy_s$(s).npz", Dict("krlist" => collect(krlist), "entropy" => ent))
end

# PLOT
plt = plot(title="Long-time average linear entropy",
    xlabel="kr (kt = kr)",
    ylabel="Linear entropy",
    legend=:bottomright)

for s in s_list
    plot!(plt, 2 .*collect(krlist), results[s], label="s = $(s)")
end

savefig(plt, "img/entropy_vs_kr_all_s.png")
display(plt)
