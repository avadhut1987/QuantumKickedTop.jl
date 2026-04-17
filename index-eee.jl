using Distributed, Plots, NPZ
addprocs()

@everywhere begin 
    using QuantumKickedTop
    using QuantumKickedTop.FloquetSystem: floquet
    using QuantumKickedTop.EigenstateEntanglement: average_entropy
end

function equipartition(N)
    if isodd(N)
        return div(N - 1, 2)
    else
        return div(N, 2)
    end
end

function eee_vs_kt(s)
    N = Int(2s)
    Q = equipartition(N)
    kr = s * π / 2

    kt_raw = LinRange(-kr, kr, 11)
    kt_set = Set(kt_raw)
    foreach(k -> push!(kt_set, k), [-kr, 0.0, kr])
    kt_vals = sort(collect(kt_set))
    avg_S_vals = pmap(kθ -> average_entropy(N, floquet(s, kr, kt); Q=Q), kt_vals)
    npzwrite("npy/eee_vs_kt_j_$(s)_partition_$(Q).npy", avg_S_vals)

    plt = plot(
        kt_vals,
        avg_S_vals,
        xlabel="kθ",
        ylabel="Normalized Von Neumann Entropy",
        title="Eigenstate Entanglement Entropy vs kθ for j = $s",
        linewidth=2,
        legend=false,
        grid=true
    )
    savefig("img/eee_vs_kt_j_$(s)_partition_$(Q).png")

end

function eee_vs_kr(s; kt = 0.0)
    N = Int(2s)
    Q = equipartition(N)

    kr_vals = collect(LinRange(0.0, 3.0, 11))

    avg_S_vals = pmap(kr -> average_entropy(N, floquet(s, kr, kt); Q=Q), kr_vals)

    npzwrite("npy/eee_vs_kr_j_$(s)_partition_$(Q).npy", avg_S_vals)

    plt = plot(
        kr_vals,
        avg_S_vals,
        xlabel="kr",
        ylabel="Normalized Von Neumann Entropy",
        title="Eigenstate Entanglement Entropy vs kr for j = $s, kt = $kt",
        linewidth=2,
        legend=false,
        grid=true
    )
    savefig("img/eee_vs_kr_j_$(s)_partition_$(Q).png")
end

s = 100
eee_vs_kr(s)
