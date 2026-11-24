using Test, QuantumOptics, QuantumKickedTop
@testset "Normalized entanglement dynamics" begin
 s=1.5; U=QuantumKickedTop.FloquetSystem.floquet(s,π/2,1.0,1.0)
 r=QuantumKickedTop.EntanglementDynamics.vn_entropy_dynamics(s,U,0.9,-0.4;rdm_dim=1,nsteps=5)
 @test r.kick==collect(0:5); @test length(r.normalized_entropy)==6
 @test r.normalized_entropy[1]<1e-10; @test all(x->-1e-12<=x<=1+1e-12,r.normalized_entropy)
 @test r.max_norm_error<1e-12; @test r.Smax_bits==1.0
end
@testset "Transient heatmap" begin
 s=1.5; U=QuantumKickedTop.FloquetSystem.floquet(s,π/2,1.0,1.0); space=(θs=[0.7,1.1],ϕs=[-0.8,0.2])
 r=QuantumKickedTop.HeatmapTimeAveragedEntropy.time_averaged_vn_entropy_grid(s,U;rdm_dim=1,nsteps=6,transient=2,space=space)
 @test size(r.entropy)==(2,2); @test r.diagnostics.post_transient_samples==4
 @test r.diagnostics.successful_points==4; @test r.diagnostics.failed_points==0
 @test all(isfinite,r.entropy); @test all(x->-1e-12<=x<=1+1e-12,r.entropy)
 direct=QuantumKickedTop.EntanglementDynamics.vn_entropy_dynamics(s,U,0.7,-0.8;rdm_dim=1,nsteps=6)
 @test isapprox(r.entropy[1,1],sum(direct.normalized_entropy[4:7])/4;atol=1e-12)
end
@testset "Partition validation" begin
 @test_throws ArgumentError QuantumKickedTop.EntanglementDynamics.maximum_vn_entropy(3,0)
 @test_throws ArgumentError QuantumKickedTop.EntanglementDynamics.maximum_vn_entropy(3,3)
end
