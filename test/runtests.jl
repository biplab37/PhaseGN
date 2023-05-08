using PhaseGN, Test

M0 = rand() + 0.5
param = Parameters()
param0 = Parameters(κ=0.0)

@testset "Helper Functions" begin
    @testset "Bisection" begin
        @test PhaseGN.bisection(x -> x^2 - 2, 0, 2) ≈ √2 atol = 1e-6
    end
    @testset "Principal Value" begin
        @test PhaseGN.PrincipalValue(0.0) == 0.0
        @test PhaseGN.PrincipalValue(1e-2, 1e-3) == 100.0
        @test PhaseGN.PrincipalValue(1e-3, 1e-2) == 0.0
    end
end

include("phase_test.jl")
include("pressure_test.jl")
include("masses_test.jl")

@testset "Mean Field" begin
    @test pressure_MF(0.01, 0.0, param) < 0.01
    @test pressure_MF(0.01, 0.0, param, norm=true) < 0.01
    @test pressure_MF(3 * rand(), 0.0, param, norm=true) < 5.0
    @test abs(pressure_MF(1.0, 0.0, param, norm=true) - pressure_MF(1.2, 0.0, param, norm=true)) < 0.1
    @test 0.0 < energy_MF(0.01, 0.0, param) < 0.01
    @test energy_MF(rand(), rand(), param) >= 0.0
    @test number_MF(rand(), rand(), param) >= 0.0
end

@testset "Parameter" begin
    @test typeof(default_parameters()) == Parameters
end
