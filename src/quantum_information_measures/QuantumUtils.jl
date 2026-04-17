module QuantumUtils

    using QuantumOptics, QuantumInterface, LinearAlgebra, SpecialFunctions, Statistics
    using ..PhiStates: get_phi_states
    export coherentspinstate, compute_rhoQ, compute_von_neumann_entropy, compute_linear_entropy

    function coherentspinstate(b, θ, φ)
        result = Ket(b)
        data = result.data
        N = length(b) - 1
        α = sin(θ / 2) * exp(1im * φ / 2)
        β = cos(θ / 2) * exp(-1im * φ / 2)

        coeff = 1.0
        factor = β^N
        @inbounds for n in 1:N+1
            data[n] = coeff * factor
            coeff *= α * sqrt((N + 1 - n) / n)          # α contribution
            factor /= β                           # progressively reduce β^power
        end
        return result
    end

    function compute_rhoQ(ψ::Vector{ComplexF64}, N::Int, Q::Int)
        A_Q = zeros(ComplexF64, Q + 1, N - Q + 1)
        log_fact = [lgamma(k + 1) for k in 0:N+1]   # real log factorial

        @inbounds for m in 0:Q
            log_comb_Q_m = log_fact[Q+1] - log_fact[m+1] - log_fact[Q-m+1]
            for n in 0:(N-Q)
                if m + n > N
                    continue
                end
                log_comb_NQ_n = log_fact[N-Q+1] - log_fact[n+1] - log_fact[N-Q-n+1]
                log_comb_N_mn = log_fact[N+1] - log_fact[m+n+1] - log_fact[N-m-n+1]
                coeff = exp(0.5 * (log_comb_Q_m + log_comb_NQ_n - log_comb_N_mn))
                A_Q[m+1, n+1] = coeff * ψ[m+n+1]
            end
        end
        return A_Q * A_Q'
    end

    function compute_von_neumann_entropy(ρ_Q::AbstractMatrix; eigval_tol::Float64=1e-12)
        ρ_Q = ComplexF64.(ρ_Q)
        trval = tr(ρ_Q)

        if isapprox(trval, 0; atol=1e-14)
            return 0.0   # entropy of the zero state = 0
        end

        ρ_Q ./= trval
        eigvals = eigen(Hermitian(ρ_Q)).values
        eigvals = eigvals[eigvals.>eigval_tol]
        return -sum(eigvals .* log2.(eigvals))
    end

    function compute_linear_entropy(ρ_Q::AbstractMatrix; eigval_tol::Float64=1e-12)
        ρ_Q = ComplexF64.(ρ_Q)
        trval = tr(ρ_Q)
        if isapprox(trval, 0; atol=1e-14)
            return 0.0
        end
        ρ_Q ./= trval
        return 1 - real(tr(ρ_Q * ρ_Q))
    end

end



# module QuantumUtils

#     using QuantumOptics, QuantumInterface, LinearAlgebra, SpecialFunctions, Statistics
#     using ..PhiStates: get_phi_states
#     export coherentspinstate, compute_rhoQ, compute_von_neumann_entropy, partial_trace_fast


#     function coherentspinstate(b, theta, phi)
#         result = Ket(b)
#         data = result.data
#         N = length(b) - 1
#         α = sin(theta / 2) * exp(1im * phi / 2)
#         β = cos(theta / 2) * exp(-1im * phi / 2)

#         coeff = 1.0
#         for n = 1:N+1
#             data[n] = coeff
#             coeff *= α * sqrt((N + 1 - n) / n)
#         end

#         factor = 1.0
#         for n = N+1:-1:1
#             data[n] *= factor
#             factor *= β
#         end
#         return result
#     end

#     function compute_rhoQ(psi::Vector{ComplexF64}, N::Int, Q::Int)
#         A_Q = zeros(ComplexF64, Q + 1, N - Q + 1)
#         log_fact = [loggamma(k + 1)[1] for k in 0:N+1]

#         for m in 0:Q
#             for n in 0:(N-Q)
#                 if m + n > N
#                     continue
#                 end
#                 log_comb_Q_m = log_fact[Q+1] - log_fact[m+1] - log_fact[Q-m+1]
#                 log_comb_NQ_n = log_fact[N-Q+1] - log_fact[n+1] - log_fact[N-Q-n+1]
#                 log_comb_N_mn = log_fact[N+1] - log_fact[m+n+1] - log_fact[N-m-n+1]
#                 log_coeff = log_comb_Q_m + log_comb_NQ_n - log_comb_N_mn
#                 coeff = exp(0.5 * log_coeff)
#                 A_Q[m+1, n+1] = coeff * psi[m+n+1]
#             end
#         end

#         return A_Q * A_Q'
#     end

#     function compute_von_neumann_entropy(rho_Q::Matrix{ComplexF64}; eigval_tol::Float64=1e-12)
#         rho_Q ./= tr(rho_Q) 
#         eigvals = eigen(Hermitian(rho_Q)).values
#         eigvals = eigvals[eigvals.>eigval_tol]
#         return -sum(eigvals .* log2.(eigvals))
#     end

#     function partial_trace_fast(rho::AbstractMatrix, traced_out::Vector{Int}, N::Int)
#         keep = setdiff(1:N, traced_out)
#         dims = fill(2, N)  # 2-level system per qubit

#         # Reshape ρ into (dims..., dims...) shape
#         ρ_tensor = reshape(rho, dims..., dims...)

#         # Axes to trace: for each traced qubit q, sum over axis q and axis N+q
#         for q in sort(traced_out; rev=true)
#             ρ_tensor = sum(ρ_tensor, dims=(q, N + q))
#         end

#         # Reshape back into matrix form
#         dim_keep = 2^(length(keep))
#         return reshape(ρ_tensor, dim_keep, dim_keep)
#     end
# end


# function rho2qubit(s::Union{Int,Float64,Rational}, ψt::AbstractVector, N::Int)
#     b = SpinBasis(s)
#     Jz = dense(0.5 * sigmaz(b))
#     Jp = 0.5 * sigmap(b)
#     vplus = real((N^2 - 2N + 4 * real(ψt' * Jz * Jz * ψt) + 4 * real(ψt' * Jz * ψt) * (N - 1)) / (4N * (N - 1)))
#     vminus = real((N^2 - 2N + 4 * real(ψt' * Jz * Jz * ψt) - 4 * real(ψt' * Jz * ψt) * (N - 1)) / (4N * (N - 1)))
#     xplus = real(((N - 1) * real(ψt' * Jp * ψt) + real(ψt' * (Jp * Jz + Jz * Jp) * ψt)) / (2N * (N - 1)))
#     xminus = real(((N - 1) * real(ψt' * Jp * ψt) - real(ψt' * (Jp * Jz + Jz * Jp) * ψt)) / (2N * (N - 1)))
#     w = real((N^2 - 4 * real(ψt' * Jz * Jz * ψt)) / (4N * (N - 1)))
#     y = w
#     u = real(ψt' * Jp * Jp * ψt) / (N * (N - 1))

#     mat = [vplus      xplus   xplus   u;
#            xplus        w       y     xminus;
#            xplus        y       w     xminus;
#            u          xminus xminus   vminus]

#     return mat
# end


# function conditional_entropy_A_measurement(psi::Ket, b_A::Basis, N::Int, subsystem_A::Int, θ::Float64, ϕ::Float64)

#     rho = compute_rhoQ(psi.data, N, N)
#     Q_B = 2^(N - subsystem_A)
#     if Q_B <= 0
#         error("Invalid subsystem split: Q_B = N - Q_A = $N - $subsystem_A = $Q_B ≤ 0")
#     end

#     b_B = NLevelBasis(2^Q_B)
#     full_basis = Operator(b_A ⊗ b_B, rho)

#     display("rho dim: ", dims(rho))
#     display("basis of the sub system A: ", b_A)
#     display("basis of the sub system B: ",  b_B)

#     ψ₊ = cos(θ / 2) * basisstate(b_A, 1) + exp(1im * ϕ) * sin(θ / 2) * basisstate(b_A, 2)
#     ψ₋ = -sin(θ / 2) * basisstate(b_A, 1) + exp(1im * ϕ) * cos(θ / 2) * basisstate(b_A, 2)

#     Π₊ = ψ₊ ⊗ dagger(ψ₊)
#     Π₋ = ψ₋ ⊗ dagger(ψ₋)

#     P₊ = embed(full_basis, subsystem_A, Π₊)
#     P₋ = embed(full_basis, subsystem_A, Π₋)

#     cond₊ = P₊ * ρ * P₊
#     cond₋ = P₋ * ρ * P₋

#     p₊ = real(tr(cond₊))
#     p₋ = real(tr(cond₋))

#     S = 0.0
#     if p₊ > 1e-12
#         ρB₊ = ptrace(cond₊ / p₊, subsystem_A)
#         S += p₊ * entropy_vn(ρB₊)
#     end
#     if p₋ > 1e-12
#         ρB₋ = ptrace(cond₋ / p₋, subsystem_A)
#         S += p₋ * entropy_vn(ρB₋)
#     end

#     return S
# end

# function min_conditional_entropy(psi::Ket, b_A::Basis, N::Int, subsystem_A::Int; θ_res::Int=20, ϕ_res::Int=20)
#     min_S = Inf
#     for θ in range(0, π, length=θ_res)
#         for ϕ in range(0, 2π, length=ϕ_res)
#             S = conditional_entropy_A_measurement(psi, b_A, N, subsystem_A, θ, ϕ)
#             if S < min_S
#                 min_S = S
#             end
#         end
#     end
#     return min_S
# end

# function entropy_matrix(ρ::Matrix{ComplexF64})
#     vals = eigen(Hermitian(ρ)).values
#     vals = clamp.(vals, 0, 1)
#     vals = vals[vals.>1e-12]
#     return -sum(vals .* log2.(vals))
# end

# function compute_discord_point(ψ::Ket, U::Operator, itr::Int, N::Int, Q_A::Int)

#     data = Array{Float64}(undef, itr)
#     for i in 1:itr
#         rho_A = compute_rhoQ(ψ.data, N, Q_A)
#         S_A = entropy_matrix(rho_A)

#         rho_B = compute_rhoQ(ψ.data, N, N - Q_A)
#         S_B = entropy_matrix(rho_B)

#         ρt = ψ ⊗ dagger(ψ)
#         S_AB = entropy_vn(ρt)

#         b_A = NLevelBasis(Q_A)
#         Φ = get_phi_states(N)
#         rho = Φ' * ρt * Φ
#         if i == 1
#             println("length of state: ", length(ψ))
#             println("rho: ", rho)
#             println("length of rho_A: ", length(rho_A))
#             println("length of rho_B: ", length(rho_B))
#         end
#         # S_cond = min_conditional_entropy(ψ, b_A, N, Q_A)

#         # data[i] = S_A + S_B - S_AB + S_cond
#         # ψ = normalize(U * ψ)
#     end
#     return mean(data)
# end