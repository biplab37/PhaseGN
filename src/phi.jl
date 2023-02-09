"""
	imagpart_ϕ(temp, μ, ω,param)

Returns the imaginary part of the polarisation for ϕ₁ and ϕ₂ at zero external momentum.
"""
function imagpart_ϕ(temp, μ, ω,param)
	m = σ1(temp,μ,param)
	Nσ = ω^2 - 4*m^2
	if Nσ > 0 && abs(ω) < 2*sqrt(param.Λ^2 + m^2)
		return ω*(1 - numberF(temp,-μ,ω/2) - numberF(temp,μ,ω/2))/(4)
	else
		return 0.0
	end
end
function imagpart_ϕ(temp, μ, ω,m,param)
	Nσ = ω^2 - 4*m^2
	if Nσ > 0 && abs(ω) < 2*sqrt(param.Λ^2 + m^2)
		return ω*(1 - numberF(temp,-μ,ω/2) - numberF(temp,μ,ω/2))/(4)
	else
		return 0.0
	end
end

"""
	realpart_ϕ(temp, μ, ω,param)

Returns the real part of the polarisation for ϕ₁ and ϕ₂ at zero external momentum.
"""
function realpart_ϕ(temp,μ,ω,param)
	m = σ1(temp,μ,param)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = - p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipalValue(ω - 2*Ep(p)) - PrincipalValue(ω + 2*Ep(p)) + PrincipalValue(Ep(p)))/π
    return integrate(integrand,0.0,param.Λ)
end
function realpart_ϕ(temp,μ,ω,m,param)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = - p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipalValue(ω - 2*Ep(p)) - PrincipalValue(ω + 2*Ep(p)) + PrincipalValue(Ep(p)))/π
    return integrate(integrand,0.0,param.Λ)
end

function fullrealpart_ϕ(temp,μ,ω,param)
	m = σ1(temp,μ,param)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = 1/π + p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipalValue(ω - 2*Ep(p)) - PrincipalValue(ω + 2*Ep(p)))/π
    return  -1/π + integrate(integrand,0.0,param.Λ)
end
function fullrealpart_ϕ(temp,μ,ω,m,param)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = 1/π + p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipalValue(ω - 2*Ep(p)) - PrincipalValue(ω + 2*Ep(p)))/π
    return  -1/π + integrate(integrand,0.0,param.Λ)
end

"""
	Π0_ϕ(temp,μ,param)

Returns the momentum and frequency independent part of the polarisation function
"""
function Π0_ϕ(temp,μ,param)
	m = σ1(temp,μ,param)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = 1/π - p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))/(π*Ep(p))
	return -1/π + integrate(integrand,0.0,param.Λ)
end

@doc raw"""
	phasesc_ϕ(temp,μ,ω,param)

Returns the scattering part of the phase for ϕ₁ and ϕ₂ at zero external momentum.

Also consider the equivalent function 

	phasesc_ϕ(temp,μ,ω,m,param)

where you supply the values of $\bar{\sigma}_1 = m$. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phasesc_ϕ(temp,μ,ω,param)
	repi = realpart_ϕ(temp,μ,ω,param)
	impi = imagpart_ϕ(temp,μ,ω,param)
	return -angle(Complex(repi,impi))
end
function phasesc_ϕ(temp,μ,ω,m,param)
	repi = realpart_ϕ(temp,μ,ω,m,param)
	impi = imagpart_ϕ(temp,μ,ω,m,param)
	return -angle(Complex(repi,impi))
end

@doc raw"""
	phaser_ϕ(temp,μ,ω,param)

Returns the resonant part of the phase for ϕ₁ and ϕ₂ at zero external momentum.

Also consider the equivalent function 

	phaser_ϕ(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency 
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases
at different vaules of frequencies at a fixed temp and μ.
"""
function phaser_ϕ(temp,μ,ω,param)
    repi = realpart_ϕ(temp,μ,ω,param)
	impi = imagpart_ϕ(temp,μ,ω,param)
	Π00 = Π0_ϕ(temp,μ,param)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end
function phaser_ϕ(temp,μ,ω,m,Π00,param)
    repi = realpart_ϕ(temp,μ,ω,m,param)
	impi = imagpart_ϕ(temp,μ,ω,m,param)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

@doc raw"""
	phase_ϕ(temp,μ,ω,param)

Returns all the phases in an array 
	[scattered phase, resonant phase, total phase]

Also consider the equivalent function 

	phase_ϕ(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phase_ϕ(temp,μ,ω,param)
    repi = realpart_ϕ(temp,μ,ω,param)
	impi = imagpart_ϕ(temp,μ,ω,param)
	Π00 = Π0_ϕ(temp,μ,param)

	phasesc = -angle(Complex(repi,impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

function phase_ϕ(temp,μ,ω,m,Π00,param)
    repi = realpart_ϕ(temp,μ,ω,m,param)
	impi = imagpart_ϕ(temp,μ,ω,m,param)

	phasesc = -angle(Complex(repi,impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

export imagpart_ϕ, realpart_ϕ, fullrealpart_ϕ, phasesc_ϕ, phaser_ϕ, phase_ϕ, Π0_ϕ
