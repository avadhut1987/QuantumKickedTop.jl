module HeatmapDiscord

    using QuantumOptics, SharedArrays, Distributed, Plots, NPZ
    import ..QuantumUtils: coherentspinstate
    import ..QuantumUtilsDiscord: compute_discord_point
    export discord_grid, density_plot_discord

    function discord_grid(s, U, Q_A, itr, space)
        b = SpinBasis(s)
        θs = space.θs
        ϕs = space.ϕs
        res_θ = length(θs)
        res_ϕ = length(ϕs)
        discord_data = SharedArray{Float64}((res_θ, res_ϕ))

        @sync @distributed for i in 1:res_θ
            for j in 1:res_ϕ
                total = 0.0
                ψt = coherentspinstate(b, θs[i], ϕs[j])
                for _ in 1:itr
                    ψt.data .= normalize(U * ψt).data
                    # compute discord for Q_A qubits
                    total += compute_discord_point(ψt.data, Q_A)
                end
                discord_data[i, j] = total / itr
            end
        end
        return discord_data
    end

    function density_plot_discord(data, Q_A, space, s, params)
        θs = space.θs
        ϕs = space.ϕs
        kr = params.kr
        kt = params.kt

        npzwrite("npy/discord_$(s)_partition_$(Q_A)_kr_$(kr)_kt_$(kt).npy", data)

        heatmap(
            collect(ϕs), collect(θs), data,
            levels=500, color=:turbo,
            xlabel="ϕ", ylabel="θ",
            title="Discord: j = $(s), Q = $(Q_A), kr = $(kr), kt = $(kt)"
        )
        savefig("img/discord_$(s)_partition_$(Q_A)_kr_$(kr)_kt_$(kt).png")
    end

end
