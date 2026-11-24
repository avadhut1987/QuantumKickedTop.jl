using Test
using QuantumOptics
using QuantumKickedTop

@testset "Normalized entanglement dynamics" begin
    s = 1.5
    Q = 1
    U = QuantumKickedTop.FloquetSystem.floquet(s, π/2, 1.0, 1.0)
    result = QuantumKickedTop.EntanglementDynamics.vn_entropy_dynamics(
        s, U, 0.9, -0.4; rdm_dim=Q, nsteps=5
    )
    @test result.kick == collect(0:5)
    @test length(result.normalized_entropy) == 6
    @test result.normalized_entropy[1] < 1e-10
    @test all(x -> -1e-12 <= x <= 1 + 1e-12, result.normalized_entropy)
    @test result.max_norm_error < 1e-12
    @test result.Smax_bits == 1.0
end

@testset "Transient-discarded normalized entropy grid" begin
    s = 1.5
    U = QuantumKickedTop.FloquetSystem.floquet(s, π/2, 1.0, 1.0)
    space = (θs=[0.7, 1.1], ϕs=[-0.8, 0.2])
    result = QuantumKickedTop.HeatmapTimeAveragedEntropy.time_averaged_vn_entropy_grid(
        s, U; rdm_dim=1, nsteps=6, transient=2, space=space
    )
    @test size(result.entropy) == (2, 2)
    @test result.diagnostics.post_transient_samples == 4
    @test result.diagnostics.successful_points == 4
    @test result.diagnostics.failed_points == 0
    @test all(isfinite, result.entropy)
    @test all(x -> -1e-12 <= x <= 1 + 1e-12, result.entropy)
end

@testset "Entropy normalization validation" begin
    @test_throws ArgumentError QuantumKickedTop.EntanglementDynamics.maximum_vn_entropy(3, 0)
    @test_throws ArgumentError QuantumKickedTop.EntanglementDynamics.maximum_vn_entropy(3, 3)
end
