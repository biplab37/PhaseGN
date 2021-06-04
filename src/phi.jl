"""
	imagpart_ϕ(temp, μ, ω,κ)

Returns the imaginary part of the polarisation for ϕ₁ and ϕ₂ at zero external momentum.
"""
function imagpart_ϕ(temp, μ, ω,κ)
	m = σ1(temp,μ,κ)
	Nσ = ω^2 - 4*m^2
	if Nσ > 0 && abs(ω) < 2*sqrt(Λ[1]^2 + m^2)
		return ω*(1 - numberF(temp,-μ,ω/2) - numberF(temp,μ,ω/2))/(4)
	else
		return 0.0
	end
end
function imagpart_ϕ(temp, μ, ω,κ,m)
	Nσ = ω^2 - 4*m^2
	if Nσ > 0 && abs(ω) < 2*sqrt(Λ[1]^2 + m^2)
		return ω*(1 - numberF(temp,-μ,ω/2) - numberF(temp,μ,ω/2))/(4)
	else
		return 0.0
	end
end

"""
	realpart_ϕ(temp, μ, ω,κ)

Returns the real part of the polarisation for ϕ₁ and ϕ₂ at zero external momentum.
"""
function realpart_ϕ(temp,μ,ω,κ)
	m = σ1(temp,μ,κ)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = - p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipleValue(ω - 2*Ep(p)) - PrincipleValue(ω + 2*Ep(p)) + PrincipleValue(Ep(p)))/π
    return hquadrature(integrand,0.0,Λ[1],reltol=1e-3,maxevals=10000)[1]
end
function realpart_ϕ(temp,μ,ω,κ,m)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = - p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipleValue(ω - 2*Ep(p)) - PrincipleValue(ω + 2*Ep(p)) + PrincipleValue(Ep(p)))/π
    return hquadrature(integrand,0.0,Λ[1],reltol=1e-3,maxevals=10000)[1]
end

function fullrealpart_ϕ(temp,μ,ω,κ)
	m = σ1(temp,μ,κ)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = 1/π + p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipleValue(ω - 2*Ep(p)) - PrincipleValue(ω + 2*Ep(p)))/π
    return  -1/π + hquadrature(integrand,0.0,Λ[1],reltol=1e-3,maxevals=10000)[1]
end
function fullrealpart_ϕ(temp,μ,ω,κ,m)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = 1/π + p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipleValue(ω - 2*Ep(p)) - PrincipleValue(ω + 2*Ep(p)))/π
    return  -1/π + hquadrature(integrand,0.0,Λ[1],reltol=1e-3,maxevals=10000)[1]
end

"""
	Π0_ϕ(temp,μ,κ)

Returns the momentum and frequency independent part of the polarisation function
"""
function Π0_ϕ(temp,μ,κ)
	m = σ1(temp,μ,κ)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = 1/π - p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))/(π*Ep(p))
	return -1/π + hquadrature(integrand,0.0,Λ[1],reltol=1e-3,maxevals=10000)[1]
end

@doc raw"""
	phasesc_ϕ(temp,μ,ω,κ)

Returns the scattering part of the phase for ϕ₁ and ϕ₂ at zero external momentum.

Also consider the equivalent function 

	phasesc_ϕ(temp,μ,ω,κ,m)

where you supply the values of $\bar{\sigma}_1 = m$. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phasesc_ϕ(temp,μ,ω,κ)
	repi = realpart_ϕ(temp,μ,ω,κ)
	impi = imagpart_ϕ(temp,μ,ω,κ)
	return angle(Complex(repi,-impi))
end
function phasesc_ϕ(temp,μ,ω,κ,m)
	repi = realpart_ϕ(temp,μ,ω,κ,m)
	impi = imagpart_ϕ(temp,μ,ω,κ,m)
	return angle(Complex(repi,-impi))
end

@doc raw"""
	phaser_ϕ(temp,μ,ω,κ)

Returns the resonant part of the phase for ϕ₁ and ϕ₂ at zero external momentum.

Also consider the equivalent function 

	phaser_ϕ(temp,μ,ω,κ,m,Π00)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency 
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases
at different vaules of frequencies at a fixed temp and μ.
"""
function phaser_ϕ(temp,μ,ω,κ)
    repi = realpart_ϕ(temp,μ,ω,κ)
	impi = imagpart_ϕ(temp,μ,ω,κ)
	Π00 = Π0_ϕ(temp,μ,κ)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end
function phaser_ϕ(temp,μ,ω,κ,m,Π00)
    repi = realpart_ϕ(temp,μ,ω,κ,m)
	impi = imagpart_ϕ(temp,μ,ω,κ,m)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

@doc raw"""
	phase_ϕ(temp,μ,ω,κ)

Returns all the phases in an array 
	[scattered phase, resonant phase, total phase]

Also consider the equivalent function 

	phase_ϕ(temp,μ,ω,κ,m,Π00)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phase_ϕ(temp,μ,ω,κ)
    repi = realpart_ϕ(temp,μ,ω,κ)
	impi = imagpart_ϕ(temp,μ,ω,κ)
	Π00 = Π0_ϕ(temp,μ,κ)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

function phase_ϕ(temp,μ,ω,κ,m,Π00)
    repi = realpart_ϕ(temp,μ,ω,κ,m)
	impi = imagpart_ϕ(temp,μ,ω,κ,m)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

export imagpart_ϕ, realpart_ϕ, fullrealpart_ϕ, phasesc_ϕ, phaser_ϕ, phase_ϕ, Π0_ϕ
