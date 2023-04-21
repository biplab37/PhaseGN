@testset "σ1" begin
    @test σ1(0.01,0.0) > 1.0 
    @test σ1(0.01,0.0,param0) ≈ 1.0 atol=1e-4
    @test σ1(0.1,0.0,param0)  < 1.0
    @test σ1(0.1,0.0) > σ1(0.1,0.1,param0)
    @test σ1(0.01,0.0,Parameters(M=M0,κ = 0.0)) ≈ M0 atol=1e-4
end

@testset "M_σ" begin
    @test M_sigma(rand(),rand(),param) >= 0.0 # mass is positive
end

@testset "M_ϕ" begin
    @test M_phi(rand(),rand(),param) >= 0.0 # mass is positive
    @test M_phi(0.1,0.0,param) <= 2.0*σ1(0.1,0.0,param) # pseudo scalar condensate is bound.
end