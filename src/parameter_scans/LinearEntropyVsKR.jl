module LinearEntropyVsKR

using QuantumOptics
using LinearAlgebra
using ..QuantumUtils: coherentspinstate, compute_rhoQ, compute_linear_entropy
using ..FloquetSystem: floquet

export long_time_avg_entropy_vs_kr

function long_time_avg_entropy_vs_kr(
    s::Float64,
    p::Float64,
    krlist::AbstractVector{<:Real},
    rdm_dim::Int,
    itr::Int;
    theta0::Real=π / 2,
    phi0::Real=-π / 2
)
    N = Int(2s)
    b = SpinBasis(s)
    ψ0 = coherentspinstate(b, theta0, phi0)

    robust_norm(x) =
        try
            sqrt(sum(abs2, x.data))
        catch
            sqrt(sum(abs2, x))
        end

    entropies = zeros(length(krlist))

    for (idx, kr) in enumerate(krlist)
        kt = kr
        k = kr + kt
        kp = kr - kt
        U = floquet(s, p, k, kp)
        ψt = deepcopy(ψ0)
        total = 0.0

        for _ = 1:itr
            ψt = U * ψt
            n = robust_norm(ψt)
            ψt.data ./= n
            ρ = compute_rhoQ(ψt.data, N, rdm_dim)
            total += compute_linear_entropy(ρ)
        end

        entropies[idx] = total / itr
        @info "s=$s, kr=$kr → entropy=$(entropies[idx])"
    end

    return entropies
end

end # module
