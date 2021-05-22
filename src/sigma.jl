function imagpart_σ(temp, μ, ω)
    m = σ1(temp,μ)
    Nσ = ω^2 - 4*m^2
    if Nσ > 0 && abs(ω) < 2*sqrt(Λ[1]^2 + m^2)
        return Nσ*(1 - numberF(temp,-μ,ω/2) - numberF(temp,μ,ω/2))/(4*ω)
    else
        return 0.0
    end
end

function Π0_σ(temp,μ)
    m = σ1(temp,μ)
    Ep(p) = sqrt(p^2+m^2)
    integrand(p) = 1/π - p*(p^2/Ep(p)^2)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))/(π*Ep(p))
    return -1/π + hquadrature(integrand,0.0,Λ[1],reltol=1e-3,maxevals=10000)[1]
end

function realpart_σ(temp,μ,ω)
    m = σ1(temp,μ)
    Ep(p) = sqrt(p^2+m^2)
    integrand(p) = - p*(p^2/Ep(p)^2)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipleValue(ω - 2*Ep(p)) - PrincipleValue(ω + 2*Ep(p)) + PrincipleValue(Ep(p)))/π
    return hquadrature(integrand,0.0,Λ[1],reltol=1e-3,maxevals=10000)[1]
end

function fullrealpart_σ(temp,μ,ω)
    m = σ1(temp,μ)
    Ep(p) = sqrt(p^2+m^2)
    integrand(p) = 1/π + p*(p^2/Ep(p)^2)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipleValue(ω - 2*Ep(p)) - PrincipleValue(ω + 2*Ep(p)))/π
    return  -1/π +hquadrature(integrand,0.0,Λ[1],reltol=1e-3,maxevals=10000)[1]
end

function phasesc_σ(temp,μ,ω,κ)
    repi = realpart_σ(temp,μ,ω,κ)
    impi = imagpart_σ(temp,μ,ω,κ)
    return angle(Complex(repi,-impi))
end

function phaser_σ(temp,μ,ω,κ)
    repi = realpart_σ(temp,μ,ω,κ)
    impi = imagpart_σ(temp,μ,ω,κ)
    Π00 = Π0_σ(temp,μ,κ)

    return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

function phase_σ(temp,μ,ω,κ)
    repi = realpart_σ(temp,μ,ω,κ)
    impi = imagpart_σ(temp,μ,ω,κ)
    Π00 = Π0_σ(temp,μ,κ)

    phasesc = angle(Complex(repi,-impi))
    phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

    return [phasesc, phaser, phasesc + phaser]
end

export imagpart_σ, realpart_σ, fullrealpart_σ, phasesc_σ, phaser_σ, phase_σ