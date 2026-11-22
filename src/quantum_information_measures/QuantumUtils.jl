
module QuantumUtils

    using QuantumOptics, QuantumInterface, LinearAlgebra, SpecialFunctions, Statistics
    using ..PhiStates: get_phi_states

    export coherentspinstate, compute_rhoQ,
           compute_von_neumann_entropy, compute_linear_entropy

    # ---------------------------------------------------------
    # Spin coherent state
    # ---------------------------------------------------------

    function coherentspinstate(b, θ, φ)

        N = length(b) - 1
        N >= 0 || throw(ArgumentError("Invalid basis dimension"))

        isfinite(θ) && isfinite(φ) ||
            throw(ArgumentError("Angles must be finite"))

        result = Ket(b)
        data = result.data

        α = sin(θ / 2) * exp(1im * φ / 2)
        β = cos(θ / 2) * exp(-1im * φ / 2)

        # Construct coefficients from the numerically
        # more stable end of the binomial distribution.

        if abs(β) >= abs(α)

            data[1] = β^N

            for k in 0:(N-1)
                data[k+2] = data[k+1] *
                    (α / β) * sqrt((N-k)/(k+1))
            end

        else

            data[N+1] = α^N

            for k in (N-1):-1:0
                data[k+1] = data[k+2] *
                    (β / α) * sqrt((k+1)/(N-k))
            end

        end

        nrm = norm(data)

        isfinite(nrm) && nrm > 0 ||
            throw(ArgumentError("Invalid coherent-state coefficients"))

        data ./= nrm

        return result
    end


    # ---------------------------------------------------------
    # Reduced density matrix for symmetric N-qubit states
    # ---------------------------------------------------------

    function compute_rhoQ(
        ψ::AbstractVector,
        N::Int,
        Q::Int
    )

        N >= 0 || throw(ArgumentError("N must be non-negative"))

        0 <= Q <= N ||
            throw(ArgumentError("Q must satisfy 0 ≤ Q ≤ N"))

        length(ψ) == N + 1 ||
            throw(ArgumentError("State vector must have length N+1"))

        all(isfinite, ψ) ||
            throw(ArgumentError("State vector contains non-finite values"))

        A_Q = zeros(ComplexF64, Q+1, N-Q+1)

        log_fact = [loggamma(k+1) for k in 0:N+1]

        @inbounds for m in 0:Q

            log_comb_Q_m =
                log_fact[Q+1] -
                log_fact[m+1] -
                log_fact[Q-m+1]

            for n in 0:(N-Q)

                log_comb_NQ_n =
                    log_fact[N-Q+1] -
                    log_fact[n+1] -
                    log_fact[N-Q-n+1]

                log_comb_N_mn =
                    log_fact[N+1] -
                    log_fact[m+n+1] -
                    log_fact[N-m-n+1]

                coeff = exp(
                    0.5 * (
                        log_comb_Q_m +
                        log_comb_NQ_n -
                        log_comb_N_mn
                    )
                )

                A_Q[m+1,n+1] = coeff * ψ[m+n+1]
            end
        end

        return A_Q * A_Q'
    end


    # ---------------------------------------------------------
    # Density matrix validation and normalization
    # ---------------------------------------------------------

    function _validated_density_matrix(
        ρ_Q::AbstractMatrix;
        herm_tol::Float64 = 1e-10,
        trace_tol::Float64 = 1e-14
    )

        size(ρ_Q,1) == size(ρ_Q,2) ||
            throw(ArgumentError("Density matrix must be square"))

        ρ = Matrix{ComplexF64}(ρ_Q)

        all(isfinite, ρ) ||
            throw(ArgumentError("Density matrix contains non-finite values"))

        trval = tr(ρ)

        abs(imag(trval)) <= herm_tol * max(abs(real(trval)), 1.0) ||
            throw(ArgumentError("Density matrix has a complex trace"))

        real(trval) > trace_tol ||
            throw(ArgumentError("Density matrix has zero or invalid trace"))

        ρ ./= real(trval)

        herm_error =
            norm(ρ - ρ') / max(norm(ρ), eps(Float64))

        herm_error <= herm_tol ||
            throw(ArgumentError("Density matrix is not Hermitian"))

        # Remove insignificant numerical asymmetry.
        ρ = (ρ + ρ') / 2

        return ρ
    end


    # ---------------------------------------------------------
    # Von Neumann entropy
    # ---------------------------------------------------------

    function compute_von_neumann_entropy(
        ρ_Q::AbstractMatrix;
        eigval_tol::Float64 = 1e-12,
        herm_tol::Float64 = 1e-10
    )

        eigval_tol >= 0 ||
            throw(ArgumentError("eigval_tol must be non-negative"))

        ρ = _validated_density_matrix(
            ρ_Q;
            herm_tol = herm_tol
        )

        λ = eigen(Hermitian(ρ)).values

        minimum(λ) >= -eigval_tol ||
            throw(ArgumentError(
                "Density matrix has significantly negative eigenvalues"
            ))

        # Clip only small negative roundoff errors.
        λ = max.(λ, 0.0)

        λ ./= sum(λ)

        # Include every positive eigenvalue, however small.
        λ = filter(x -> x > 0.0, λ)

        return -sum(λ .* log2.(λ))
    end


    # ---------------------------------------------------------
    # Linear entropy
    # ---------------------------------------------------------

    function compute_linear_entropy(
        ρ_Q::AbstractMatrix;
        eigval_tol::Float64 = 1e-12,
        herm_tol::Float64 = 1e-10
    )

        eigval_tol >= 0 ||
            throw(ArgumentError("eigval_tol must be non-negative"))

        ρ = _validated_density_matrix(
            ρ_Q;
            herm_tol = herm_tol
        )

        λ = eigen(Hermitian(ρ)).values

        minimum(λ) >= -eigval_tol ||
            throw(ArgumentError(
                "Density matrix has significantly negative eigenvalues"
            ))

        λ = max.(λ, 0.0)
        λ ./= sum(λ)

        return 1.0 - sum(abs2, λ)
    end

end
