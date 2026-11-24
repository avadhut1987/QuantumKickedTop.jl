module HeatmapTimeAveragedEntropy

using Base.Threads
using Statistics
using ..EntanglementDynamics: vn_entropy_dynamics

export time_averaged_vn_entropy_grid

"""
    time_averaged_vn_entropy_grid(s, U; rdm_dim, nsteps, transient, space)

Compute normalized time-averaged von Neumann entropy on a coherent-state
grid. Kick zero is the product-state initial condition. The transient
kicks 0:transient are discarded, and the average uses exactly
transient+1:nsteps. Failed grid points remain NaN and are identified in
the returned diagnostics. The early/late difference is a finite-time
stability diagnostic, not evidence of asymptotic convergence.
"""
function time_averaged_vn_entropy_grid(
    s, U;
    rdm_dim::Integer,
    nsteps::Integer,
    transient::Integer,
    space
)
    nsteps >= 2 || throw(ArgumentError("nsteps must be at least 2"))
    0 <= transient < nsteps ||
        throw(ArgumentError("transient must satisfy 0 ≤ transient < nsteps"))
    θs, ϕs = space.θs, space.ϕs
    !isempty(θs) && !isempty(ϕs) ||
        throw(ArgumentError("The phase-space grid must be non-empty"))
    all(x -> x isa Real && isfinite(x) && 0 <= x <= π, θs) ||
        throw(ArgumentError("theta grid values must be finite and lie in [0, π]"))
    all(x -> x isa Real && isfinite(x) && -π <= x <= π, ϕs) ||
        throw(ArgumentError("phi grid values must be finite and lie in [-π, π]"))

    nθ, nϕ = length(θs), length(ϕs)
    average = fill(NaN, nθ, nϕ)
    early_late_difference = fill(NaN, nθ, nϕ)
    failed = falses(nθ, nϕ)
    npost = nsteps - transient
    split = max(1, fld(npost, 2))

    @threads for linear_index in 1:(nθ * nϕ)
        i = fld(linear_index - 1, nϕ) + 1
        j = mod(linear_index - 1, nϕ) + 1
        try
            result = vn_entropy_dynamics(
                s, U, θs[i], ϕs[j];
                rdm_dim=rdm_dim, nsteps=nsteps
            )
            post = @view result.normalized_entropy[(transient + 2):(nsteps + 1)]
            length(post) == npost ||
                error("Unexpected post-transient sample count")
            average[i, j] = mean(post)
            if npost >= 2
                early_late_difference[i, j] =
                    abs(mean(@view post[1:split]) - mean(@view post[(split + 1):end]))
            end
        catch err
            failed[i, j] = true
            @warn "Time-averaged entanglement failed at grid point ($i,$j)" exception=(err, catch_backtrace())
        end
    end

    diagnostics = (
        successful_points = count(!, failed),
        failed_points = count(failed),
        post_transient_samples = npost,
        early_late_difference = early_late_difference,
        failed_mask = failed
    )
    return (entropy = average, diagnostics = diagnostics)
end

end # module
