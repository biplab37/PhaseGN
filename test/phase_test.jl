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