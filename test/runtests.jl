using PhaseGN, Test

M0 = rand() + 0.5
param = Parameters()

@testset "Helper Functions" begin
    @testset "Bisection" begin
        @test PhaseGN.bisection(x->x^2-2,0,2) ≈ √2 atol=1e-6
    end
    @testset "Principal Value" begin
        @test PhaseGN.PrincipalValue(0.0) == 0.0
        @test PhaseGN.PrincipalValue(1e-2,1e-3) == 100.0
        @test PhaseGN.PrincipalValue(1e-3,1e-2) == 0.0
    end
end

@testset "Scalar Channel (Mean Field)" begin
    @testset "σ1" begin
        @test σ1(0.01,0.0) > 1.0 
        @test σ1(0.01,0.0,Parameters(κ=0.0)) ≈ 1.0 atol=1e-4
        @test σ1(0.1,0.0,Parameters(κ=0.0))  < 1.0
        @test σ1(0.1,0.0) > σ1(0.1,0.1,Parameters(κ = 0.0))
        @test σ1(0.01,0.0,Parameters(M=M0,κ = 0.0)) ≈ M0 atol=1e-4
    end
end

@testset "Scalar Channel" begin
    @test imagpart_sigma(0.01,0.0,0.1,param) ≈ 0.0 atol=1e-3

    @test -π <= phasesc_sigma(rand(),rand(),rand(),param) <= 0.0 
    @test 0.0 <= phaser_sigma(rand(),rand(),rand(),param) <= π
    @test phase_sigma(0.1,0.1,0.1,param)[1:2] == [phasesc_sigma(0.1,0.1,0.1,param), phaser_sigma(0.1,0.1,0.1,param)]
    @test phase_sigma(rand(),rand(),rand(),param)[3] >= 0.0
end

@testset "Pseudo-Scalar Channel" begin
    @test imagpart_phi(0.01,0.0,0.1,param) ≈ 0.0 atol=1e-3

    @test -π <= phasesc_phi(rand(),rand(),100*rand(),param) <= 0.0
    @test 0.0 <= phaser_phi(rand(),rand(),100*rand(),param) <= π
    @test phase_phi(0.1,0.1,0.1,param)[1:2] == [phasesc_phi(0.1,0.1,0.1,param), phaser_phi(0.1,0.1,0.1,param)]
    @test phase_phi(rand(),rand(),100*rand(),Parameters(κ=0.1))[3] >= 0.0
end

@testset "Pressure" begin
    include("pressure_test.jl")
end

@testset "Mean Field" begin
    @test pressure_MF(0.01,0.0,param) < 0.001
    @test pressure_MF(0.01,0.0,param) < 0.01
    @test 0.0< energy_MF(0.01,0.0,param) < 0.01
    @test energy_MF(rand(),rand(),param) >= 0.0
    @test number_MF(rand(),rand(),param) >= 0.0 
end

@testset "Parameter" begin
    @test typeof(default_parameters()) == Parameters
end

@testset "Phase" begin
    include("phase_test.jl")
end