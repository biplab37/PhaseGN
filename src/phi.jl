"""
	imagpart_phi(temp, μ, ω,param)

Returns the imaginary part of the polarisation for ϕ₁ and ϕ₂ at zero external momentum.
"""
function imagpart_phi(ω, temp, μ, param)
    m = σ1(temp, μ, param)
    Nσ = ω^2 - 4 * m^2
    if Nσ > 0 && abs(ω) < 2 * sqrt(param.Λ^2 + m^2)
        return ω * (1 - numberF(temp, -μ, ω / 2) - numberF(temp, μ, ω / 2)) / (4)
    else
        return 0.0
    end
end

function imagpart_phi(ω, temp, μ, m, param)
    Nσ = ω^2 - 4 * m^2
    if Nσ > 0 && abs(ω) < 2 * sqrt(param.Λ^2 + m^2)
        return ω * (1 - numberF(temp, -μ, ω / 2) - numberF(temp, μ, ω / 2)) / (4)
    else
        return 0.0
    end
end

"""
	realpart_phi(temp, μ, ω,param)

Returns the real part of the polarisation for ϕ₁ and ϕ₂ at zero external momentum.
"""
function realpart_phi(ω, temp, μ, param)
    m = σ1(temp, μ, param)
    Ep(p) = sqrt(p^2 + m^2)
    integrand(p) = -p * (1) * (1 - numberF(temp, -μ, Ep(p)) - numberF(temp, μ, Ep(p))) * (PrincipalValue(ω - 2 * Ep(p)) - PrincipalValue(ω + 2 * Ep(p)) + PrincipalValue(Ep(p))) / π
    return integrate(integrand, 0.0, param.Λ)
end

function realpart_phi(ω, temp, μ, m, param)
    Ep(p) = sqrt(p^2 + m^2)
    integrand(p) = -p * (1) * (1 - numberF(temp, -μ, Ep(p)) - numberF(temp, μ, Ep(p))) * (PrincipalValue(ω - 2 * Ep(p)) - PrincipalValue(ω + 2 * Ep(p)) + PrincipalValue(Ep(p))) / π
    return integrate(integrand, 0.0, param.Λ)
end

function fullrealpart_phi(ω, temp, μ, param)
    m = σ1(temp, μ, param)
    Ep(p) = sqrt(p^2 + m^2)
    integrand(p) = 1 / π + p * (1) * (1 - numberF(temp, -μ, Ep(p)) - numberF(temp, μ, Ep(p))) * (PrincipalValue(ω - 2 * Ep(p)) - PrincipalValue(ω + 2 * Ep(p))) / π
    return -1 / π + integrate(integrand, 0.0, param.Λ)
end

function fullrealpart_phi(ω, temp, μ, m, param)
    Ep(p) = sqrt(p^2 + m^2)
    integrand(p) = 1 / π + p * (1) * (1 - numberF(temp, -μ, Ep(p)) - numberF(temp, μ, Ep(p))) * (PrincipalValue(ω - 2 * Ep(p)) - PrincipalValue(ω + 2 * Ep(p))) / π
    return -1 / π + integrate(integrand, 0.0, param.Λ)
end

"""
	Π0_phi(temp,μ,param)

Returns the momentum and frequency independent part of the polarisation function
"""
function Π0_phi(temp, μ, param)
    m = σ1(temp, μ, param)
    Ep(p) = sqrt(p^2 + m^2)
    integrand(p) = 1 / π - p * (1) * (1 - numberF(temp, -μ, Ep(p)) - numberF(temp, μ, Ep(p))) / (π * Ep(p))
    return -1 / π + integrate(integrand, 0.0, param.Λ)
end

function Π0_phi_m(temp, μ, m, param)
    Ep(p) = sqrt(p^2 + m^2)
    integrand(p) = 1 / π - p * (1) * (1 - numberF(temp, -μ, Ep(p)) - numberF(temp, μ, Ep(p))) / (π * Ep(p))
    return -1 / π + integrate(integrand, 0.0, param.Λ)
end

@doc raw"""
	phasesc_phi(temp,μ,ω,param)

Returns the scattering part of the phase for ϕ₁ and ϕ₂ at zero external momentum.

Also consider the equivalent function 

	phasesc_phi(temp,μ,ω,m,param)

where you supply the values of $\bar{\sigma}_1 = m$. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phasesc_phi(ω, temp, μ, param)
    repi = realpart_phi(ω, temp, μ, param)
    impi = imagpart_phi(ω, temp, μ, param)
    return -angle(Complex(repi, impi))
end

function phasesc_phi(ω, temp, μ, m, param)
    repi = realpart_phi(ω, temp, μ, m, param)
    impi = imagpart_phi(ω, temp, μ, m, param)
    return -angle(Complex(repi, impi))
end

@doc raw"""
	phaser_phi(temp,μ,ω,param)

Returns the resonant part of the phase for ϕ₁ and ϕ₂ at zero external momentum.

Also consider the equivalent function 

	phaser_phi(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency 
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases
at different vaules of frequencies at a fixed temp and μ.
"""
function phaser_phi(ω, temp, μ, param)
    repi = realpart_phi(ω, temp, μ, param)
    impi = imagpart_phi(ω, temp, μ, param)
    Π00 = Π0_phi(temp, μ, param)

    return -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))
end
function phaser_phi(ω, temp, μ, m, Π00, param)
    repi = realpart_phi(ω, temp, μ, m, param)
    impi = imagpart_phi(ω, temp, μ, m, param)

    return -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))
end

@doc raw"""
	phase_phi(ω,temp,μ,param)

Returns all the phases in an array 
	[scattered phase, resonant phase, total phase]

Also consider the equivalent function 

	phase_phi(ω,temp,μ,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phase_phi(ω, temp, μ, param)
    repi = realpart_phi(ω, temp, μ, param)
    impi = imagpart_phi(ω, temp, μ, param)
    Π00 = Π0_phi(temp, μ, param)

    phasesc = -angle(Complex(repi, impi))
    phaser = -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))

    return [phasesc, phaser, phasesc + phaser]
end

function phase_phi(ω, temp, μ, m, Π00, param)
    repi = realpart_phi(ω, temp, μ, m, param)
    impi = imagpart_phi(ω, temp, μ, m, param)

    phasesc = -angle(Complex(repi, impi))
    phaser = -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))

    return [phasesc, phaser, phasesc + phaser]
end

function mass_phi(temp, μ, param)
    func(ME) = fullrealpart_phi(ME, temp, μ, param)
    return bisection(func, 0.0, param.Λ)
end

export imagpart_phi, realpart_phi, fullrealpart_phi, phasesc_phi, phaser_phi, phase_phi, Π0_phi
