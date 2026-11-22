module HeatmapEntropy

using QuantumOptics, SharedArrays, Distributed, Plots, NPZ, LinearAlgebra
using ..QuantumUtils: coherentspinstate, compute_rhoQ, compute_von_neumann_entropy
using ..FloquetSystem: floquet

export vn_entropy_grid, density_plot

function _validate_floquet_operator(U, basis; tolerance=1e-10)

    hasproperty(U, :data) ||
        throw(ArgumentError("U must be a QuantumOptics operator with matrix data"))

    hasproperty(U, :basis_l) && hasproperty(U, :basis_r) ||
        throw(ArgumentError("U must expose left and right Hilbert-space bases"))

    U.basis_l == basis ||
        throw(ArgumentError("Floquet operator left basis does not match SpinBasis(s)"))

    U.basis_r == basis ||
        throw(ArgumentError("Floquet operator right basis does not match SpinBasis(s)"))

    M = U.data
    D = length(basis)

    M isa AbstractMatrix ||
        throw(ArgumentError("Floquet operator data must be a matrix"))

    size(M) == (D, D) ||
        throw(DimensionMismatch(
            "Floquet operator must have size ($D, $D); got $(size(M))"
        ))

    all(isfinite, M) ||
        throw(ArgumentError("Floquet operator contains non-finite entries"))

    residual = norm(adjoint(M) * M - I) / sqrt(D)

    isfinite(residual) && residual <= tolerance ||
        throw(ArgumentError(
            "Floquet operator is not unitary within tolerance. " *
            "Normalized unitarity residual = $residual"
        ))

    return nothing
end


function vn_entropy_grid(
    s, U::Any, rdm_dim, itr, space;
    recompute_U_on_worker::Bool=false,
    floquet_args=nothing,
    max_bad_per_point::Int=5
)
    # ---------------- INPUT VALIDATION ----------------

    s isa Real && isfinite(s) && s >= 0 ||
        throw(ArgumentError("s must be a finite, non-negative real number"))

    isapprox(2s, round(2s); atol=1e-12, rtol=0) ||
        throw(ArgumentError("s must be an integer or half-integer"))

    N = Int(round(2s))
    rdm_dim isa Integer || throw(ArgumentError("rdm_dim must be an integer"))
    0 <= rdm_dim <= N || throw(ArgumentError("rdm_dim must satisfy 0 ≤ rdm_dim ≤ N"))
    itr isa Integer && itr > 0 ||
        throw(ArgumentError("itr must be a positive integer"))
    max_bad_per_point >= 0 ||
        throw(ArgumentError("max_bad_per_point must be non-negative"))

    recompute_U_on_worker && floquet_args === nothing &&
        throw(ArgumentError("recompute_U_on_worker=true but no floquet_args provided"))

    b = SpinBasis(s)

    # Validate the Floquet operator before starting the grid calculation.
    # In worker-reconstruction mode, construct one copy here for preflight.
    U_check = recompute_U_on_worker ? floquet(floquet_args...) : U
    _validate_floquet_operator(U_check, b)

    θs, ϕs = space.θs, space.ϕs

    res_θ, res_ϕ = length(θs), length(ϕs)

    res_θ > 0 && res_ϕ > 0 ||
        throw(ArgumentError("The phase-space grid must be non-empty"))

    all(x -> x isa Real && isfinite(x), θs) ||
        throw(ArgumentError("θ grid must contain only finite real values"))

    all(x -> x isa Real && isfinite(x), ϕs) ||
        throw(ArgumentError("ϕ grid must contain only finite real values"))

    all(θ -> 0 <= θ <= π, θs) ||
        throw(ArgumentError("θ values must lie in [0, π]"))

    all(ϕ -> -π <= ϕ <= π, ϕs) ||
        throw(ArgumentError("ϕ values must lie in [-π, π]"))

    # Initialize all shared arrays explicitly.
    # NaN means that the entropy has not been successfully computed.
    entropy = SharedArray{Float64}((res_θ, res_ϕ))
    diagnostics = SharedArray{UInt8}((res_θ, res_ϕ))

    fill!(entropy, NaN)
    fill!(diagnostics, 0)

    # ---------------- SINGLE-TRAJECTORY CALCULATION ----------------

    function calculate_point(U_local, θ, ϕ)
        ψt = coherentspinstate(b, θ, ϕ)
        total = 0.0

        for kick in 1:itr
            # A failed propagation invalidates the complete trajectory.
            ψt_temp = U_local * ψt

            n = norm(ψt_temp.data)

            isfinite(n) && n > 0.0 ||
                error("Invalid state norm at kick $kick: $n")

            # Normalize the propagated state.
            ψt.data .= ψt_temp.data ./ n

            ρQ = compute_rhoQ(ψt.data, N, Int(rdm_dim))
            ent = compute_von_neumann_entropy(ρQ)

            isfinite(ent) ||
                error("Non-finite entropy at kick $kick")

            total += ent
        end

        average_entropy = total / itr

        isfinite(average_entropy) ||
            error("Non-finite time-averaged entropy")

        return average_entropy
    end

    # ---------------- FIRST PARALLEL PASS ----------------
    #
    # Assign rows explicitly to workers. Each active worker constructs
    # its Floquet operator once and reuses it for all assigned rows.

    worker_ids = workers()

    isempty(worker_ids) && (worker_ids = [myid()])

    n_active = min(length(worker_ids), res_θ)
    active_worker_ids = worker_ids[1:n_active]

    # Round-robin row assignment balances the number of rows per worker.
    row_chunks = [
        collect(worker_index:n_active:res_θ)
        for worker_index in 1:n_active
    ]

    @sync for (worker_index, pid) in enumerate(active_worker_ids)

        assigned_rows = row_chunks[worker_index]

        @async remotecall_wait(pid, assigned_rows) do rows

            # Construct or receive the Floquet operator only once
            # for this worker's complete assignment.
            U_local = try
                if recompute_U_on_worker
                    U_worker = floquet(floquet_args...)
                    _validate_floquet_operator(U_worker, b)
                    U_worker
                else
                    U
                end
            catch err
                @warn "Could not initialize or validate Floquet operator on worker $pid" exception=(err, catch_backtrace())
                nothing
            end

            for i in rows
                for j in 1:res_ϕ

                    if U_local === nothing
                        entropy[i, j] = NaN
                        diagnostics[i, j] = 1
                        continue
                    end

                    try
                        entropy[i, j] = calculate_point(
                            U_local, θs[i], ϕs[j]
                        )
                        diagnostics[i, j] = 0

                    catch err
                        # Keep failed trajectories distinguishable from zero entropy.
                        entropy[i, j] = NaN
                        diagnostics[i, j] = 1

                        @warn "Entropy trajectory failed at grid point ($i,$j)" exception=(err, catch_backtrace())
                    end
                end
            end
        end
    end

    # ---------------- SECOND PASS: RETRY FAILED POINTS ----------------

    n_bad = count(x -> x == 1, Array(diagnostics))

    println("First pass completed.")
    println("Points requiring recomputation = $n_bad")

    if n_bad > 0 && max_bad_per_point > 0

        bad_points = [
            (i, j)
            for i in 1:res_θ, j in 1:res_ϕ
            if diagnostics[i, j] == 1
        ]

        println("Retrying $(length(bad_points)) failed points...")

        # Reuse the operator already constructed and validated during preflight.
        U_retry = U_check

        if U_retry !== nothing

            for (i, j) in bad_points

                success = false
                last_error = nothing

                for retry in 1:max_bad_per_point

                    try
                        # Each retry starts from the original coherent state.
                        value = calculate_point(
                            U_retry, θs[i], ϕs[j]
                        )

                        entropy[i, j] = value
                        diagnostics[i, j] = 0
                        success = true

                        println("Recovered grid point ($i,$j) on retry $retry")
                        break

                    catch err
                        last_error = err
                    end
                end

                if !success
                    # NaN is retained: failure is not zero entropy.
                    entropy[i, j] = NaN
                    diagnostics[i, j] = 1

                    @warn "Grid point ($i,$j) failed after $max_bad_per_point retries; retaining NaN" exception=last_error
                end
            end
        end
    end

    # ---------------- FINAL DIAGNOSTICS ----------------

    n_bad_after = count(x -> x == 1, Array(diagnostics))

    println("Entropy calculation finished.")
    println("Successfully computed points = $(res_θ * res_ϕ - n_bad_after)")
    println("Unresolved failed points = $n_bad_after")

    return entropy, "points_with_failed_trajectories = $n_bad_after"
end


function density_plot(data, rdm_dim, space, s, params)
    θs, ϕs = space.θs, space.ϕs
    p, k, kp = params.p, params.k, params.kp

    npzwrite(
        "npy/entropy_j$(s)_partition$(rdm_dim)_p$(p)_k$(k)_kp$(kp).npy",
        data
    )

    heatmap(
        collect(ϕs), collect(θs), data,
        levels=500,
        color=:turbo,
        xlabel="ϕ",
        ylabel="θ",
        title="VN: j = $(s), k = $(k), kp = $(kp)"
    )

    savefig(
        "img/entropy_j$(s)_partition$(rdm_dim)_p$(p)_k$(k)_kp$(kp).png"
    )
end

end # module
