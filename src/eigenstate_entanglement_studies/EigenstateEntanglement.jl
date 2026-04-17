module EigenstateEntanglement
    using QuantumOptics, LinearAlgebra, Statistics, SharedArrays, Distributed, Plots, NPZ
    import ..QuantumUtils: coherentspinstate, compute_rhoQ, compute_von_neumann_entropy
    export average_entropy

    
    function average_entropy(N::Int, U::Operator; Q=div(N - 1, 2))

        s = N / 2
        _, eigenvecs = eigen(U.data)

        vn_entropies = map(i -> begin
                psi = eigenvecs[:, i]
                rho_Q = compute_rhoQ(psi, N, Q)
                compute_von_neumann_entropy(rho_Q)
            end, 1:(N+1))

        return mean(vn_entropies) / log2(Q + 1)
    end
end