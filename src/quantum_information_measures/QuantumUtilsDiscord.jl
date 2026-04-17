module QuantumUtilsDiscord

    using LinearAlgebra
    using ..QuantumUtils: compute_von_neumann_entropy, compute_rhoQ
    export extract_subsystem, compute_discord_rho, compute_discord_point, reduce_density_matrix

    function extract_subsystem(ψ::AbstractVector, N::Int, keep::Vector{Int})
        @assert length(keep) == 1 "compute_rhoQ only supports one index"
        return compute_rhoQ(ψ, N, keep[1])
    end

    function partial_trace_fast(ρ::AbstractMatrix, traced_out::Vector{Int}, N::Int)
        dims = ntuple(_ -> 2, N)
        ρ_tensor = reshape(ρ, dims..., dims...)

        # Little-endian: qubit 1 is least significant
        for q in sort(traced_out; rev=true)
            ρ_tensor = sum(ρ_tensor, dims=(N - q + 1, 2N - q + 1))
        end

        keep = setdiff(1:N, traced_out)
        dim_keep = 2^length(keep)

        return reshape(ρ_tensor, dim_keep, dim_keep)
    end

    function reduce_density_matrix(ρ::AbstractMatrix, N::Int, keep::Vector{Int})
        traced_out = setdiff(1:N, keep)
        return partial_trace_fast(ρ, traced_out, N)
    end

    function computational_projectors(Q_A::Int)
        dim = 2^Q_A
        return [Diagonal([i == j ? 1.0 : 0.0 for i in 1:dim]) for j in 1:dim]
    end


    function compute_discord_rho(ρ_AB::AbstractMatrix, Q_A::Int; dθ=π / 30, dφ=π / 30)
        Ns = Int(round(log2(size(ρ_AB, 1))))  # total number of qubits
        Q_B = Ns - Q_A
        idx_A = collect(1:Q_A)
        idx_B = collect(Q_A+1:Ns)

        # --- Step 1: Reduced density matrices ---
        ρ_A = reduce_density_matrix(ρ_AB, Ns, idx_A)
        ρ_B = reduce_density_matrix(ρ_AB, Ns, idx_B)

        # --- Step 2: Entropies ---
        S_A = compute_von_neumann_entropy(ρ_A)
        S_B = compute_von_neumann_entropy(ρ_B)
        S_AB = compute_von_neumann_entropy(ρ_AB)

        # --- Step 3: Classical correlations ---
        min_cond_entropy = Inf

        if Q_A == 1
            I_B = I(2^Q_B)

            θ_vals = 0:dθ:π
            φ_vals = 0:dφ:(2π-dφ)

            for θ in θ_vals, φ in φ_vals
                c, s = cos(θ / 2), sin(θ / 2)
                phase = exp(im * φ)

                ket0 = [c, phase * s]
                ket1 = [-s, phase * c]
                Π0 = ket0 * ket0'
                Π1 = ket1 * ket1'

                cond_entropy = 0.0
                for Π in (Π0, Π1)
                    M = kron(Π, I_B)
                    ρ_proj = M * ρ_AB * M
                    p_k = real(tr(ρ_proj))
                    if p_k > 1e-12
                        ρ_Bk = reduce_density_matrix(ρ_proj, Ns, idx_B) / p_k
                        cond_entropy += p_k * compute_von_neumann_entropy(ρ_Bk)
                    end
                end
                min_cond_entropy = min(min_cond_entropy, cond_entropy)
            end
        else
            # fallback: computational basis projectors
            basis = computational_projectors(Q_A)
            I_B = I(2^Q_B)
            cond_entropy = 0.0
            for Π in basis
                M = kron(Π, I_B)
                ρ_proj = M * ρ_AB * M
                p_k = real(tr(ρ_proj))
                if p_k > 1e-12
                    ρ_Bk = reduce_density_matrix(ρ_proj, Ns, idx_B) / p_k
                    cond_entropy += p_k * compute_von_neumann_entropy(ρ_Bk)
                end
            end
            min_cond_entropy = cond_entropy
        end

        # --- Step 4: Quantum Discord ---
        C_A = S_B - min_cond_entropy
        I_AB = S_A + S_B - S_AB
        return I_AB - C_A
    end

    function compute_discord_point(ψ::AbstractVector, Q_A::Int; dθ=π / 30, dφ=π / 30)
        ρ_AB = ψ * ψ'
        return compute_discord_rho(ρ_AB, Q_A; dθ=dθ, dφ=dφ)
    end

    function compute_discord_point(ψ::AbstractVector, N::Int, Ns::Int, Q_A::Int; dθ=π / 30, dφ=π / 30)
        keep = collect(1:Ns)
        ρ_sub = reduce_density_matrix(ψ * ψ', N, keep)
        return compute_discord_rho(ρ_sub, Q_A; dθ=dθ, dφ=dφ)
    end

end
