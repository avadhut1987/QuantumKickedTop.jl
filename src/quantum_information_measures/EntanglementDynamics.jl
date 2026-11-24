module EntanglementDynamics
using QuantumOptics, LinearAlgebra
using ..QuantumUtils: coherentspinstate, compute_rhoQ, compute_von_neumann_entropy
export vn_entropy_dynamics, maximum_vn_entropy, validate_floquet_operator
function maximum_vn_entropy(N::Integer, Q::Integer)
    N >= 2 || throw(ArgumentError("N must be at least 2"))
    1 <= Q < N || throw(ArgumentError("rdm_dim must satisfy 1 ≤ Q < N"))
    log2(min(Q + 1, N - Q + 1))
end
function validate_floquet_operator(U, basis; tolerance=1e-10)
    hasproperty(U, :data) && hasproperty(U, :basis_l) && hasproperty(U, :basis_r) || throw(ArgumentError("U must be a QuantumOptics operator"))
    U.basis_l == basis && U.basis_r == basis || throw(ArgumentError("Floquet operator basis mismatch"))
    D = length(basis)
    size(U.data) == (D,D) || throw(DimensionMismatch("Floquet operator must have size ($D, $D)"))
    all(isfinite, U.data) || throw(ArgumentError("Floquet operator contains non-finite entries"))
    residual = norm(adjoint(U.data)*U.data-I)/sqrt(D)
    isfinite(residual) && residual <= tolerance || throw(ArgumentError("Floquet operator is not unitary; residual=$residual"))
    nothing
end
"""
    vn_entropy_dynamics(s, U, theta0, phi0; rdm_dim, nsteps)
Return normalized von Neumann entanglement entropy at kicks 0:nsteps.
The coherent product state at kick zero is included. Normalization is
log2(min(Q+1,N-Q+1)) bits for the Q | (N-Q) symmetric partition.
"""
function vn_entropy_dynamics(s, U, theta0::Real, phi0::Real; rdm_dim::Integer, nsteps::Integer, validate_operator::Bool=true)
    s isa Real && isfinite(s) && s >= 0 || throw(ArgumentError("s must be finite and non-negative"))
    isapprox(2s,round(2s);atol=1e-12,rtol=0) || throw(ArgumentError("s must be integer or half-integer"))
    N=Int(round(2s)); Q=Int(rdm_dim); Smax=maximum_vn_entropy(N,Q)
    nsteps >= 0 || throw(ArgumentError("nsteps must be non-negative"))
    isfinite(theta0) && 0 <= theta0 <= π || throw(ArgumentError("theta0 must lie in [0, π]"))
    isfinite(phi0) && -π <= phi0 <= π || throw(ArgumentError("phi0 must lie in [-π, π]"))
    basis=SpinBasis(s); validate_operator && validate_floquet_operator(U,basis)
    ψ=coherentspinstate(basis,theta0,phi0)
    entropy=Vector{Float64}(undef,nsteps+1); state_norm=similar(entropy)
    for kick in 0:nsteps
        nrm=norm(ψ.data); isfinite(nrm) && nrm>0 || error("Invalid state norm at kick $kick")
        state_norm[kick+1]=nrm; ψ.data ./= nrm
        ρA=compute_rhoQ(ψ.data,N,Q); value=compute_von_neumann_entropy(ρA)
        isfinite(value) || error("Non-finite entropy at kick $kick")
        entropy[kick+1]=clamp(value/Smax,0.0,1.0)
        kick<nsteps && (ψ=U*ψ)
    end
    (kick=collect(0:nsteps),normalized_entropy=entropy,state_norm=state_norm,Smax_bits=Smax,max_norm_error=maximum(abs.(state_norm.-1.0)))
end
end
