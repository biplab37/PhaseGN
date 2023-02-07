using PhaseGN, Test

M0 = rand()
param = Parameters()

@testset "Helper Functions" begin
    @testset "Bisection" begin
        @test PhaseGN.bisection(x->x^2-2,0,2) ≈ √2 atol=1e-6
    end
    @testset "Principal Value" begin
        @test PhaseGN.PrincipalValue(1e-2,1e-3) == 100.0
        @test PhaseGN.PrincipalValue(1e-3,1e-2) == 0.0
    end
end

@testset "Scalar Channel (Mean Field)" begin
    @testset "σ1" begin
        @test σ1(0.01,0.0) ≈ 1.0 atol=1e-4
        @test σ1(0.01,0.0,Parameters(κ=0.01)) > 1.0
        @test σ1(0.1,0.0)  < 1.0
        @test σ1(0.1,0.0) > σ1(0.1,0.1,Parameters())
        @test σ1(0.01,0.0,Parameters(M=M0)) ≈ M0 atol=1e-4
    end
end

@testset "Scalar Channel" begin
    @test imagpart_σ(0.01,0.0,0.1,param) ≈ 0.0 atol=1e-3

    @test -π <= phasesc_σ(rand(),rand(),rand(),param) <= π 
    @test -π <= phaser_σ(rand(),rand(),rand(),param) <= π
    @test phase_σ(0.1,0.1,0.1,param)[1:2] == [phasesc_σ(0.1,0.1,0.1,param), phaser_σ(0.1,0.1,0.1,param)]
end

@testset "Pseudo-Scalar Channel" begin
    @test imagpart_ϕ(0.01,0.0,0.1,param) ≈ 0.0 atol=1e-3

    @test -π <= phasesc_ϕ(rand(),rand(),rand(),param) <= π 
    @test -π <= phaser_ϕ(rand(),rand(),rand(),param) <= π
    @test phase_ϕ(0.1,0.1,0.1,param)[1:2] == [phasesc_ϕ(0.1,0.1,0.1,param), phaser_ϕ(0.1,0.1,0.1,param)]   
end

## Future Works, currently gettting skipped

@testset "Mean Field" begin
    @test pressure_MF(0.01,0.0,param) < 0.01 skip=true
    @test energy_MF(0.01,0.0,param) < 0.01 skip=true
    @test number_MF(rand(),rand(),param) > 0.0 skip=true
end