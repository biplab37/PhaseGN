temp = rand()
μ = 0.0
param = Parameters()

@testset "External Momentum" begin
    @testset "Pseudo-Scalar Channel" begin
        ph_sc, ph_r, ph_tot = phase_phi(temp, μ, 10 * rand(), 10 * rand(), param)
        @test -π <= ph_sc <= 0.0
        @test 0.0 <= ph_r <= π
    end
    @testset "Scalar Channel" begin
        si_sc, si_r, si_tot = phase_sigma(temp, μ, 10 * rand(), 10 * rand(), param)
        @test -π <= si_sc <= 0.0
        @test 0.0 <= si_r <= π
    end
end
