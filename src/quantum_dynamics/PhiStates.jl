module PhiStates

    using Combinatorics, LinearAlgebra
    export get_phi_states

    # Build a W-like state with `q` excitations
    function W_state(n::Int, q::Int)
        basis = collect(permutations(vcat(zeros(Int, n - q), ones(Int, q))))
        norm = sqrt(binomial(n, q))
        state = sum([kron([b == 0 ? [1, 0] : [0, 1] for b in bvec]...) for bvec in basis]) / norm
        return state
    end

    function Wbar_state(n::Int, q::Int)
        basis = collect(permutations(vcat(ones(Int, q), zeros(Int, n - q))))
        norm = sqrt(binomial(n, q))
        state = sum([kron([b == 1 ? [1, 0] : [0, 1] for b in bvec]...) for bvec in basis]) / norm
        return state
    end

    # Φ⁺ state
    function phi_plus(n::Int, q::Int)
        phase = im^(n - 2q)
        (W_state(n, q) + phase * Wbar_state(n, q)) / sqrt(2)
    end

    # Φ⁻ state
    function phi_minus(n::Int, q::Int)
        phase = im^(n - 2q)
        (W_state(n, q) - phase * Wbar_state(n, q)) / sqrt(2)
    end

    # Exported function to get all Φ states
    function get_phi_states(n::Int)
        Φ = []

        for q in 0:div(n, 2)
            if isodd(n) || q < div(n, 2)
                push!(Φ, phi_plus(n, q))
            else
                push!(Φ, W_state(n, q))
            end
        end

        if isodd(n)
            for q in 0:div(n - 1, 2)
                push!(Φ, phi_minus(n, q))
            end
        else
            for q in 0:div(n, 2) - 1
                push!(Φ, phi_minus(n, q))
            end
        end

        return Φ
    end

end # module
