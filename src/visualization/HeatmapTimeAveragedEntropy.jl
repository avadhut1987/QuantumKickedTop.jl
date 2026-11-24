module HeatmapTimeAveragedEntropy
using Base.Threads, Statistics, QuantumOptics
using ..EntanglementDynamics: vn_entropy_dynamics, validate_floquet_operator
export time_averaged_vn_entropy_grid
"""Compute normalized entropy averages over kicks transient+1:nsteps.
Kicks 0:transient are discarded. Failed points remain NaN. Early/late
differences are finite-time stability diagnostics, not convergence proofs."""
function time_averaged_vn_entropy_grid(s,U;rdm_dim::Integer,nsteps::Integer,transient::Integer,space)
 nsteps>=2 || throw(ArgumentError("nsteps must be at least 2"))
 0<=transient<nsteps || throw(ArgumentError("transient must satisfy 0 ≤ transient < nsteps"))
 θs,ϕs=space.θs,space.ϕs
 !isempty(θs)&&!isempty(ϕs) || throw(ArgumentError("Grid must be non-empty"))
 all(x->x isa Real&&isfinite(x)&&0<=x<=π,θs) || throw(ArgumentError("theta values must lie in [0, π]"))
 all(x->x isa Real&&isfinite(x)&&-π<=x<=π,ϕs) || throw(ArgumentError("phi values must lie in [-π, π]"))
 validate_floquet_operator(U,SpinBasis(s))
 nθ,nϕ=length(θs),length(ϕs); avg=fill(NaN,nθ,nϕ); delta=fill(NaN,nθ,nϕ); failed=falses(nθ,nϕ)
 npost=nsteps-transient; split=max(1,fld(npost,2))
 @threads for idx in 1:(nθ*nϕ)
  i=fld(idx-1,nϕ)+1; j=mod(idx-1,nϕ)+1
  try
   r=vn_entropy_dynamics(s,U,θs[i],ϕs[j];rdm_dim=rdm_dim,nsteps=nsteps,validate_operator=false)
   post=@view r.normalized_entropy[(transient+2):(nsteps+1)]
   length(post)==npost || error("Unexpected post-transient sample count")
   avg[i,j]=mean(post)
   npost>=2 && (delta[i,j]=abs(mean(@view(post[1:split]))-mean(@view(post[(split+1):end]))))
  catch err
   failed[i,j]=true
   @warn "Entropy failed at grid point ($i,$j)" exception=(err,catch_backtrace())
  end
 end
 (entropy=avg,diagnostics=(successful_points=count(!,failed),failed_points=count(failed),post_transient_samples=npost,early_late_difference=delta,failed_mask=failed))
end
end
