module ClassicalUtils

using LinearAlgebra, Statistics, StaticArrays
using ..ClassicalMap: map, tangent_map

export spherical_to_cartesian, get_LLE, get_lyapunov_spectrum, ks_entropy

@inline function spherical_to_cartesian(theta, phi)
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)
    return SVector{3,Float64}(x, y, z)
end

# Construct an orthonormal basis for the tangent plane at R.
@inline function tangent_basis(R::SVector{3,Float64})
    reference = abs(R[3]) < 0.9 ?
        SVector(0.0, 0.0, 1.0) :
        SVector(1.0, 0.0, 0.0)

    e1 = normalize(cross(R, reference))
    e2 = normalize(cross(R, e1))

    return hcat(e1, e2)
end

"""
    get_lyapunov_spectrum(space, n, params; burnin=0)

Calculate the two Lyapunov exponents of the DKT classical map
using tangent-plane dynamics and periodic QR orthonormalization.

`space` must provide `theta` and `phi`.
`params` must provide `p`, `k`, and `kprime`.

Returns a two-element vector sorted in descending order.
"""
function get_lyapunov_spectrum(
    space,
    n::Int,
    params::NamedTuple;
    burnin::Int=0
)
    n > 0 || throw(ArgumentError("n must be positive"))
    0 <= burnin < n ||
        throw(ArgumentError("burnin must satisfy 0 <= burnin < n"))

    R = spherical_to_cartesian(space.theta, space.phi)

    # Orthonormal tangent basis and initial QR orientation.
    B = Matrix(tangent_basis(R))
    Q = Matrix{Float64}(I, 2, 2)

    log_sums = zeros(Float64, 2)

    for t in 1:n
        J = Matrix(tangent_map(R, params))

        Rnext = map(R, params)
        Bnext = Matrix(tangent_basis(Rnext))

        # Jacobian represented between consecutive tangent planes.
        A = transpose(Bnext) * J * B

        # Propagate and orthonormalize tangent vectors.
        V = A * Q
        F = qr(V)

        Q = Matrix(F.Q)
        triangular = Matrix(F.R)

        if t > burnin
            diagonal = abs.(diag(triangular))
            any(x -> !isfinite(x) || x <= 0.0, diagonal) &&
                throw(ErrorException("Invalid QR growth factor at step $t"))

            log_sums .+= log.(diagonal)
        end

        R = Rnext
        B = Bnext
    end

    exponents = log_sums ./ (n - burnin)
    return sort(exponents, rev=true)
end

"""
    get_LLE(space, n, params; burnin=0)

Return the largest Lyapunov exponent.
"""
function get_LLE(space, n::Int, params::NamedTuple; burnin::Int=0)
    return get_lyapunov_spectrum(space, n, params; burnin=burnin)[1]
end

"""
    ks_entropy(space, n, params; burnin=0)

Return the sum of positive Lyapunov exponents (the Pesin
estimate of the Kolmogorov-Sinai entropy per kick).
"""
function ks_entropy(space, n::Int, params::NamedTuple; burnin::Int=0)
    exponents = get_lyapunov_spectrum(space, n, params; burnin=burnin)
    return sum(max.(exponents, 0.0))
end

end
