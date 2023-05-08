@testset "Pressure" begin
    @test pressure(phasesc_sigma, 0.1, 0.0, param) <= 0.0
    @test pressure(phaser_sigma, 0.1, 0.0, param) >= 0.0
    @test pressure(phasesc_phi, 0.1, 0.0, param) <= 0.0
    @test pressure(phaser_phi, 0.1, 0.0, param) >= 0.0
end
