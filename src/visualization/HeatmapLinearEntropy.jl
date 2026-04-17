module HeatmapLinearEntropy

using QuantumOptics, SharedArrays, Distributed, Plots, NPZ, LinearAlgebra
using ..QuantumUtils: coherentspinstate, compute_rhoQ, compute_linear_entropy
export linear_entropy_grid, density_plot

function linear_entropy_grid(
    s, U::Any, rdm_dim, itr, space;
    recompute_U_on_worker::Bool=false,
    floquet_args=nothing,
    max_bad_per_point::Int=5)
    N = Int(2s)
    b = SpinBasis(s)
    θs, ϕs = space.θs, space.ϕs
    res_θ, res_ϕ = length(θs), length(ϕs)
    entropy = SharedArray{Float64}((res_θ, res_ϕ))
    diagnostics = SharedArray{UInt8}((res_θ, res_ϕ))

    function robust_norm_arr(x)
        a = try
            x.data
        catch
            x
        end
        return sqrt(sum(abs2, a))
    end

    # ---------------- FIRST PARALLEL PASS ----------------
    @sync @distributed for i in 1:res_θ
        U_local = U
        if recompute_U_on_worker
            if floquet_args === nothing
                error("recompute_U_on_worker=true but no floquet_args provided")
            end
            U_local = floquet(floquet_args...)
        end

        for j in 1:res_ϕ
            total = 0.0
            bad_count = 0
            ψt = coherentspinstate(b, θs[i], ϕs[j])

            for _ in 1:itr
                ψt_temp = try
                    U_local * ψt
                catch
                    bad_count += 1
                    continue
                end

                n = robust_norm_arr(ψt_temp)
                if isnan(n) || isinf(n) || n == 0.0
                    diagnostics[i, j] = 1
                    bad_count += 1
                    ψt = coherentspinstate(b, π / 2, 0.0)
                    continue
                end

                ψt.data .= ψt_temp.data ./ n
                ρQ = compute_rhoQ(ψt.data, N, rdm_dim)
                ent = compute_linear_entropy(ρQ)
                total += isnan(ent) || isinf(ent) ? 0.0 : ent
            end

            entropy[i, j] = total / itr
            if bad_count > 0
                diagnostics[i, j] = 1
            end
        end
    end

    # ---------------- SECOND PASS: RECOMPUTE INVALID POINTS ----------------
    n_bad = count(x -> x == 1, Array(diagnostics))
    println("Diagnostics: points_with_invalid_norms = $n_bad")

    if n_bad > 0
        bad_points = [(i, j) for i in 1:res_θ, j in 1:res_ϕ if diagnostics[i, j] == 1]
        println("Recomputing $(length(bad_points)) invalid points...")

        for (i, j) in bad_points
            retry = 0
            success = false

            while retry < max_bad_per_point && !success
                retry += 1
                try
                    total = 0.0
                    ψt = coherentspinstate(b, θs[i], ϕs[j])
                    for _ in 1:itr
                        ψt_temp = U * ψt
                        n = robust_norm_arr(ψt_temp)
                        if isnan(n) || isinf(n) || n == 0.0
                            ψt = coherentspinstate(b, π / 2, 0.0)
                            continue
                        end
                        ψt.data .= ψt_temp.data ./ n
                        ρQ = compute_rhoQ(ψt.data, N, rdm_dim)
                        ent = compute_linear_entropy(ρQ)
                        total += isnan(ent) || isinf(ent) ? 0.0 : ent
                    end
                    entropy[i, j] = total / itr
                    diagnostics[i, j] = 0
                    success = true
                catch
                    # keep retrying if it fails
                end
            end

            if !success
                @warn "Point ($i,$j) failed after $max_bad_per_point retries; keeping 0.0"
                entropy[i, j] = 0.0
            end
        end

        n_bad_after = count(x -> x == 1, Array(diagnostics))
        println("Recomputation finished. Remaining bad points = $n_bad_after")
    end

    return entropy, "points_with_invalid_norms = $(count(x -> x == 1, Array(diagnostics)))"
end


function density_plot(data, rdm_dim, space, s, params)
    θs, ϕs = space.θs, space.ϕs
    kr, kt = params.kr, params.kt
    npzwrite("npy/entropy_$(s)_partition_$(rdm_dim)_kr_$(kr)_kt_$(kt).npy", data)
    heatmap(collect(ϕs), collect(θs), data, levels=500, color=:turbo,
        xlabel="ϕ", ylabel="θ",
        title="VN: j = $(s), kr = $(kr), kt = $(kt)")
    savefig("img/entropy_$(s)_partition_$(rdm_dim)_kr_$(kr)_kt_$(kt).png")
end

end # module
