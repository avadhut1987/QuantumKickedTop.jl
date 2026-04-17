module HeatmapLLE
    using Distributed, SharedArrays, LinearAlgebra, NPZ, Plots
    using Base.Threads: Atomic, atomic_add!
    using ..ClassicalUtils: get_LLE

    export lle_grid, density_plot_lle

    function lle_grid(space, params, time, res; max_bad_per_point::Int=5)
        θs, ϕs = space.θs, space.ϕs
        lle = SharedArray{Float64}((res, res))  # Use Float64 for consistency
        diagnostics = SharedArray{UInt8}((res, res))  # 0=good, 1=bad

        @sync @distributed for idx in 1:(res*res)
            i = ((idx - 1) ÷ res) + 1
            j = ((idx - 1) % res) + 1

            total_bad = 0
            value = NaN

            # First pass: try computation with retries
            while total_bad ≤ max_bad_per_point
                try
                    value = get_LLE((theta=θs[i], phi=ϕs[j]), time, params)
                    if isfinite(value)
                        lle[i, j] = value
                        diagnostics[i, j] = 0
                        break
                    else
                        total_bad += 1
                    end
                catch
                    total_bad += 1
                end
            end

            # Mark as bad if still failed
            if !isfinite(lle[i, j])
                lle[i, j] = 0.0  # Use 0.0 instead of NaN for plotting
                diagnostics[i, j] = 1
            end
        end

        # Diagnostics summary
        n_bad = count(x -> x == 1, Array(diagnostics))
        println("Diagnostics: points_with_invalid_LLE = $n_bad")

        return lle
    end

    function density_plot_lle(data, space, params)
        θs, ϕs = space.θs, space.ϕs
        k, kprime = params.k, params.kprime

        mkpath.(["npy/", "img/"])
        npzwrite("npy/lle_k_$(round(k; digits=2))_kprime_$(round(kprime; digits=2)).npy", data)

        p = heatmap(collect(ϕs), collect(θs), data, levels=500, color=:turbo, size=(800, 600), xlabel="φ", ylabel="θ")
        savefig(p, "img/lle_k_$(round(k; digits=2))_kprime_$(round(kprime; digits=2)).png")
        println("✓ Plot saved!")
    end


end
