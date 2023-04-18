"""
	imagpart_sigma(temp, μ, ω,param)

Returns the imaginary part of the polarisation for σ₁ and σ₂ at zero external momentum.
"""
function imagpart_sigma(temp, μ, ω,param::Parameters)
    m = σ1(temp,μ,param)
    Nσ = ω^2 - 4*m^2
    if Nσ > 0 && abs(ω) < 2*sqrt(param.Λ^2 + m^2)
        return Nσ*(1 - numberF(temp,-μ,ω/2) - numberF(temp,μ,ω/2))/(4*ω)
    else
        return 0.0
    end
end

function imagpart_sigma(temp, μ, ω,m,param::Parameters)
    Nσ = ω^2 - 4*m^2
    if Nσ > 0 && abs(ω) < 2*sqrt(param.Λ^2 + m^2)
        return Nσ*(1 - numberF(temp,-μ,ω/2) - numberF(temp,μ,ω/2))/(4*ω)
    else
        return 0.0
    end
end

"""
	Π0_sigma(temp,μ,param)

Returns the momentum and frequency independent part of the polarisation function
"""
function Π0_sigma(temp,μ,param)
    m = σ1(temp,μ,param)
    Ep(p) = sqrt(p^2+m^2)
    integrand(p) = 1/π - p*(p^2/Ep(p)^2)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))/(π*Ep(p))
    return -1/π + integrate(integrand,0.0,param.Λ)
end

"""
	realpart_sigma(temp, μ, ω,param)

Returns the real part of the polarisation for σ₁ and σ₂ at zero external momentum.
"""
function realpart_sigma(temp,μ,ω,param)
    m = σ1(temp,μ,param)
    Ep(p) = sqrt(p^2+m^2)
    integrand(p) = - p*(p^2/Ep(p)^2)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipalValue(ω - 2*Ep(p)) - PrincipalValue(ω + 2*Ep(p)) + PrincipalValue(Ep(p)))/π
    return integrate(integrand,0.0,param.Λ)
end

function realpart_sigma(temp,μ,ω,m,param)
    Ep(p) = sqrt(p^2+m^2)
    integrand(p) = - p*(p^2/Ep(p)^2)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipalValue(ω - 2*Ep(p)) - PrincipalValue(ω + 2*Ep(p)) + PrincipalValue(Ep(p)))/π
    return integrate(integrand,0.0,param.Λ)
end

function fullrealpart_sigma(temp,μ,ω,param)
    m = σ1(temp,μ,param)
    Ep(p) = sqrt(p^2+m^2)
    integrand(p) = 1/π + p*(p^2/Ep(p)^2)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipalValue(ω - 2*Ep(p)) - PrincipalValue(ω + 2*Ep(p)))/π
    return  -1/π +integrate(integrand,0.0,param.Λ)
end

function fullrealpart_sigma(temp,μ,ω,m,param)
    Ep(p) = sqrt(p^2+m^2)
    integrand(p) = 1/π + p*(p^2/Ep(p)^2)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipalValue(ω - 2*Ep(p)) - PrincipalValue(ω + 2*Ep(p)))/π
    return  -1/π +integrate(integrand,0.0,param.Λ)
end

@doc raw"""
	phasesc_sigma(temp,μ,ω,param)

Returns the scattering part of the phase for σ₁ and σ₂ at zero external momentum.

Also consider the equivalent function 

    phasesc_phi(temp,μ,ω,m,param)

where you supply the values of $\bar{\sigma}_1 = m$. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phasesc_sigma(temp,μ,ω,param)
    repi = realpart_sigma(temp,μ,ω,param)
    impi = imagpart_sigma(temp,μ,ω,param)
    return -angle(Complex(repi,impi))
end

function phasesc_sigma(temp,μ,ω,m,param)
    repi = realpart_sigma(temp,μ,ω,m,param)
    impi = imagpart_sigma(temp,μ,ω,m,param)
    return -angle(Complex(repi,impi))
end

@doc raw"""
	phaser_sigma(temp,μ,ω,param)

Returns the resonant part of the phase for σ₁ and σ₂ at zero external momentum.

Also consider the equivalent function 

    phaser_phi(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency 
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases
at different vaules of frequencies at a fixed temp and μ.
"""
function phaser_sigma(temp,μ,ω,param)
    repi = realpart_sigma(temp,μ,ω,param)
    impi = imagpart_sigma(temp,μ,ω,param)
    Π00 = Π0_sigma(temp,μ,param)

    return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

function phaser_sigma(temp,μ,ω,m,Π00,param)
    repi = realpart_sigma(temp,μ,ω,m,param)
    impi = imagpart_sigma(temp,μ,ω,m,param)

    return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

@doc raw"""
	phase_sigma(temp,μ,ω,param)

Returns all the phases in an array 

	[scattered phase, resonant phase, total phase]

Also consider the equivalent function 

    phase_sigma(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phase_sigma(temp,μ,ω,param)
    repi = realpart_sigma(temp,μ,ω,param)
    impi = imagpart_sigma(temp,μ,ω,param)
    Π00 = Π0_sigma(temp,μ,param)

    phasesc = -angle(Complex(repi,impi))
    phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

    return [phasesc, phaser, phasesc + phaser]
end

function phase_sigma(temp,μ,ω,m,Π00,param)
    repi = realpart_sigma(temp,μ,ω,m,param)
    impi = imagpart_sigma(temp,μ,ω,m,param)

    phasesc = -angle(Complex(repi,impi))
    phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

    return [phasesc, phaser, phasesc + phaser]
end

export imagpart_sigma, realpart_sigma, fullrealpart_sigma, phasesc_sigma, phaser_sigma, phase_sigma,Π0_sigma
