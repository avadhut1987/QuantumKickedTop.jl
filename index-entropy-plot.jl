# index-entropy-plot.jl

using QuantumOptics
using Plots
using NPZ
using QuantumKickedTop.FloquetSystem: floquet
import QuantumKickedTop.QuantumUtils: coherentspinstate, compute_rhoQ, compute_von_neumann_entropy


function entropy_vs_time(s, U, θ0, ϕ0, nsteps, rdm_dim)
    N = Int(2s)
    b = SpinBasis(s)
    ψt = coherentspinstate(b, θ0, ϕ0)

    entropy = zeros(Float64, nsteps + 1)

    # n = 0 → entropy of initial state
    ρQ0 = compute_rhoQ(ψt.data, N, rdm_dim)
    entropy[1] = compute_von_neumann_entropy(ρQ0)

    # n = 1…nsteps → apply U
    for n in 1:nsteps
        ψt = normalize(U * ψt)
        ρQ = compute_rhoQ(ψt.data, N, rdm_dim)
        entropy[n+1] = compute_von_neumann_entropy(ρQ)
    end
    return entropy
end

function entropy_time_plot(entropy, s, params; θ0, ϕ0, rdm_dim)
    kr, kt = params.kr, params.kt
    nsteps = length(entropy) - 1

    # npzwrite("npy/vn_entropy_vs_time_j$(s)_Q$(rdm_dim)_kr$(kr)_kt$(kt).npy", entropy)

    plot(
        0:nsteps, entropy,
        xlabel="Discrete time step n",
        ylabel="Von Neumann Entropy",
        lw=2,
        title="VN Entropy vs Time: j=$(s), kr=$(kr), kt=$(kt), θ=$(θ0), ϕ=$(ϕ0), Q=$(rdm_dim)",
        legend=false,
        xticks=0:1:nsteps
    )

    savefig("img/vn_entropy_vs_time_j$(s)_Q$(rdm_dim)_kr$(kr)_kt$(kt).png")
end


function main()
    sLIst = [100.5, 200.5, 500.5]
    θ0, ϕ0 = π/2, 0.0
    nsteps = 1000
    rdm_dim = 1
    p = 2π
    k = 0.5
    kprimeList = LinRange(0.01, 6.5, 51)
    
    colors = [:red, :green, :blue]
    plt = plot(
        xlabel="k′",
        ylabel="S",
        lw=3.5,
        legend=:topright
    )

    for (i, s) in enumerate(sLIst)
        println("Running for s = $s")
        N = Int(2s)

        entropyList = Array{Float64}(undef, length(kprimeList))

        for (idx, kprime) in enumerate(kprimeList)
            U = floquet(s, p, k, kprime)
            ψt = coherentspinstate(SpinBasis(s), θ0, ϕ0)

            for _ in 1:nsteps
                ψt = normalize(U * ψt)
            end
            ρQ = compute_rhoQ(ψt.data, N, rdm_dim)
            ent = compute_von_neumann_entropy(ρQ)
            entropyList[idx] = ent
        end

        filename = "npy/entropy_vs_kprime_j$(s)_Q$(rdm_dim)_k$(k).npy"
        npzwrite(filename, entropyList)

        plot!(
            plt,
            kprimeList,
            entropyList,
            label="j = $s",
            color=colors[i]
        )
    end
    savefig(plt, "img/entropy_vs_kprime_p0.png")
end

main()




# function entropy_time_plot(entropy, s, params; θ0, ϕ0, rdm_dim)
#     kr, kt = params.kr, params.kt
#     nsteps = length(entropy) - 1

#     # npzwrite("npy/vn_entropy_vs_time_j$(s)_Q$(rdm_dim)_kr$(kr)_kt$(kt).npy", entropy)

#     plot(
#         0:nsteps, entropy,
#         xlabel="Discrete time step n",
#         ylabel="Von Neumann Entropy",
#         lw=2,
#         title="VN Entropy vs Time: j=$(s), kr=$(kr), kt=$(kt), θ=$(θ0), ϕ=$(ϕ0), Q=$(rdm_dim)",
#         legend=false,
#         xticks=0:1:nsteps
#     )

#     savefig("img/vn_entropy_vs_time_j$(s)_Q$(rdm_dim)_kr$(kr)_kt$(kt).png")
# end

# function main()
#     θ0, ϕ0 = π/2, -π/2
#     nsteps = 20
#     rdm_dim = 1

#     for s in [75.5]
#         kr = (s * π / 4)
#         kt = 100.0
#         params = (kr=kr, kt=kt)

#         println("Running for s=$s, kr=$kr, kt=$kt")

#         U = floquet(s, kr, kt)
#         entropy_data = entropy_vs_time(s, U, θ0, ϕ0, nsteps, rdm_dim)

#         entropy_time_plot(entropy_data, s, params; θ0=θ0, ϕ0=ϕ0, rdm_dim=rdm_dim)
#     end
# end

# main()