using PhaseGN, Test

M0 = rand()

@testset "Scalar Channel Mean Field" begin
    @test σ1(0.01,0.0) ≈ 1.0 atol=1e-4
    @test σ1(0.01,0.0,Parameters(κ=0.01)) > 1.0
    @test σ1(0.1,0.0)  < 1.0
    @test σ1(0.1,0.0) > σ1(0.1,0.1,Parameters())
    @test σ1(0.01,0.0,Parameters(M=M0)) ≈ M0 atol=1e-4
end