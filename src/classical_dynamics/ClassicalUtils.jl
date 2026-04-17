module ClassicalUtils

    using LinearAlgebra, Statistics, ForwardDiff, StaticArrays
    using ..ClassicalMap: map, tangent_map

    export spherical_to_cartesian, get_LLE, get_lyapunov_spectrum, ks_entropy


    @inline function spherical_to_cartesian(theta, phi)
        x = sin(theta) * cos(phi)
        y = sin(theta) * sin(phi)
        z = cos(theta)
        return SVector{3,Float64}(x, y, z)
    end

    @inline @inline function get_LLE(space, n::Int, params::NamedTuple)
        R = spherical_to_cartesian(space.theta, space.phi)
        M = Matrix{BigFloat}(I, 3, 3)
        # sumlog = log(norm(mat, 2))

        for _ in 1:n
            J = Matrix{BigFloat}(tangent_map(R, params))
            M = J * M
            R = map(R, params)
            # Z = J * M
            # M, Rmat = qr(Z)
            # sumlog += log(abs(Rmat[1,1]))
            if any(isnan.(M)) || any(isinf.(M))
                println("NaN or Inf found in mat at iteration $i")
            end
        end
        eigenvalues = eigvals(M)

        return log(sort(abs.(eigenvalues), rev=true)[1]) / n
    end

end