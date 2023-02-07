# This file contains the code for the case with non zero external momentum

#-------------------------------------------------------------#
#|                      Pseudo Scalar(ϕ)                     |#
#-------------------------------------------------------------#


"""
	imagpart_ϕ_q(temp,μ,ω,q,param)

Returns the imaginary part of the polarisation function at finite external momentum q for the pseudo scalar channel.
"""
function imagpart_ϕ_q(temp,μ,ω,q,param)
	β = 1/temp
	m = σ1(temp,μ,param)
	s = ω^2 - q^2
	Ep(p) = sqrt(p^2 + m^2)

	if s>=0
		if s<=4*m^2 || ω > 2*param.Λ
			result = 0.0
		else
			p1 = abs(-q + ω*sqrt(1 - 4*m^2/s))/2
			p2 = (q + ω*sqrt(1 - 4*m^2/s))/2
			integrand(p) = s*(numberF(temp,μ,-ω + Ep(p)) - numberF(temp,μ,Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s - 2*ω*Ep(p))/(2*p*q))^2))

			result = hquadrature(integrand,p1,p2,reltol=1e-2,maxevals=10000)[1]
		end
	else 
		p1 = abs(-q + ω*sqrt(1 - 4*m^2/s))/2
		p2 = (q + ω*sqrt(1 - 4*m^2/s))/2
		integrand1(p) = -s*(numberF(temp,μ,-ω - Ep(p)) - numberF(temp,μ,-Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s + 2*ω*Ep(p))/(2*p*q))^2))
		result = hquadrature(integrand1,p1,p2,reltol=1e-2,maxevals=10000)[1]

		integrand2(p) = -s*(numberF(temp,μ,-ω - Ep(p)) - numberF(temp,μ,-Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s + 2*ω*Ep(p))/(2*p*q))^2)) - s*(numberF(temp,μ,-ω + Ep(p)) - numberF(temp,μ,Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s - 2*ω*Ep(p))/(2*p*q))^2))
		# int1(p) = integrand2(1/(1 - p))/(1 - p)^2

		result += hquadrature(integrand2,p2,param.Λ,reltol=1e-2,maxevals=10000)[1]
	end

	return result
end
function imagpart_ϕ_q(temp,μ,ω,q,m,param)
	β = 1/temp
	s = ω^2 - q^2
	Ep(p) = sqrt(p^2 + m^2)

	if s>=0
		if s<=4*m^2 || ω > 2*param.Λ
			result = 0.0
		else
			p1 = abs(-q + ω*sqrt(1 - 4*m^2/s))/2
			p2 = (q + ω*sqrt(1 - 4*m^2/s))/2
			integrand(p) = s*(numberF(temp,μ,-ω + Ep(p)) - numberF(temp,μ,Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s - 2*ω*Ep(p))/(2*p*q))^2))

			result = hquadrature(integrand,p1,p2,reltol=1e-2,maxevals=10000)[1]
		end
	else 
		p1 = abs(-q + ω*sqrt(1 - 4*m^2/s))/2
		p2 = (q + ω*sqrt(1 - 4*m^2/s))/2
		integrand1(p) = -s*(numberF(temp,μ,-ω - Ep(p)) - numberF(temp,μ,-Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s + 2*ω*Ep(p))/(2*p*q))^2))
		result = hquadrature(integrand1,p1,p2,reltol=1e-2,maxevals=10000)[1]

		integrand2(p) = -s*(numberF(temp,μ,-ω - Ep(p)) - numberF(temp,μ,-Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s + 2*ω*Ep(p))/(2*p*q))^2)) - s*(numberF(temp,μ,-ω + Ep(p)) - numberF(temp,μ,Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s - 2*ω*Ep(p))/(2*p*q))^2))
		int1(p) = integrand2(1/(1 - p))/(1 - p)^2

		result += hquadrature(int1,1 - 1/p2,1,reltol=1e-2,maxevals=10000)[1]
	end

	return result
end


"""
	realpart_ϕ_q(temp, μ, ω, q, param)

Returns the real part of the polarisation function at finite external momentum for the pseudo scalar channel.
Uses `imagpart_ϕ_q` and Kramer-Cronig Relation for the calculation
"""
function realpart_ϕ_q(temp, μ, ω, q, param)
	
	integrand(ν) = 2*ν*imagpart_ϕ_q(temp,μ,ν,q,param)*(PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2))/π
	int1(ν) = integrand(1/(1 - ν))/(1 - ν)^2

	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]  + hquadrature(int1,0,1,reltol=1e-2,maxevals=10000)[1]
end

function realpart_ϕ_q(temp, μ, ω, q, m, param)
	
	integrand(ν) = 2*ν*imagpart_ϕ_q(temp,μ,ν,q,m,param)*(PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2))/π
	int1(ν) = integrand(1/(1 - ν))/(1 - ν)^2

	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]  + hquadrature(int1,0,1,reltol=1e-2,maxevals=10000)[1]
end

@doc raw"""
	phasesc_ϕ_q(temp,μ,ω,param)

Returns the scattering part of the phase for ϕ₁ and ϕ₂ at external momentum q.

Also consider the equivalent function 

	phasesc_ϕ_q(temp,μ,ω,m,param)

where you supply the values of $\bar{\sigma}_1 = m$. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phasesc_ϕ_q(temp,μ,ω,q,param)
	repi = realpart_ϕ_q(temp,μ,ω,q,param)
	impi = imagpart_ϕ_q(temp,μ,ω,q,param)
	return angle(Complex(repi,-impi))
end
function phasesc_ϕ_q(temp,μ,ω,q,m,param)
	repi = realpart_ϕ_q(temp,μ,ω,q,m,param)
	impi = imagpart_ϕ_q(temp,μ,ω,q,m,param)
	
	return angle(Complex(repi,-impi))
end

@doc raw"""
	phaser_ϕ_q(temp,μ,ω,param)

Returns the resonant part of the phase for ϕ₁ and ϕ₂ at external momentum q.

Also consider the equivalent function 

	phaser_ϕ_q(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency 
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases
at different vaules of frequencies at a fixed temp and μ.
"""
function phaser_ϕ_q(temp,μ,ω,q,param)
    repi = realpart_ϕ_q(temp,μ,ω,q,param)
	impi = imagpart_ϕ_q(temp,μ,ω,q,param)
	Π00 = Π0_ϕ(temp,μ,param)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end
function phaser_ϕ_q(temp,μ,ω,q,m,Π00,param)
    repi = realpart_ϕ_q(temp,μ,ω,q,m,param)
	impi = imagpart_ϕ_q(temp,μ,ω,q,m,param)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

@doc raw"""
	phase_ϕ_q(temp,μ,ω,param)

Returns all the phases for ϕ₁ and ϕ₂ at external momentum q in an array 
	[scattered phase, resonant phase, total phase]

Also consider the equivalent function 

	phase_ϕ_q(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phase_ϕ_q(temp,μ,ω,q,param)
    repi = realpart_ϕ_q(temp,μ,ω,q,param)
	impi = imagpart_ϕ_q(temp,μ,ω,q,param)
	Π00 = Π0_ϕ(temp,μ,param)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

function phase_ϕ_q(temp,μ,ω,q,m,Π00,param)
    repi = realpart_ϕ_q(temp,μ,ω,q,m,param)
	impi = imagpart_ϕ_q(temp,μ,ω,q,m,param)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

#-------------------------------------------------------------#
#|                         Scalar(σ)                         |#
#-------------------------------------------------------------#

"""
	imagpart_σ_q(temp, μ, ω,param)

Returns the imaginary part of the polarisation for σ₁ and σ₂ at external momentum q.
"""
function imagpart_σ_q(temp,μ,ω,q,param)
	β = 1/temp
	m = σ1(temp,μ,param)
	s = ω^2 - q^2
	Ep(p) = sqrt(p^2 + m^2)

	if s>=0
		if s<=4*m^2 || ω > 2*param.Λ
			result = 0.0
		else
			p1 = abs(-q + ω*sqrt(1 - 4*m^2/s))/2
			p2 = (q + ω*sqrt(1 - 4*m^2/s))/2
			integrand(p) = (s-4*m^2)*(numberF(temp,μ,-ω + Ep(p)) - numberF(temp,μ,Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s - 2*ω*Ep(p))/(2*p*q))^2))

			result = hquadrature(integrand,p1,p2,reltol=1e-2,maxevals=10000)[1]
		end
	else 
		p1 = abs(-q + ω*sqrt(1 - 4*m^2/s))/2
		p2 = (q + ω*sqrt(1 - 4*m^2/s))/2
		integrand1(p) = -(s-4*m^2)*(numberF(temp,μ,-ω - Ep(p)) - numberF(temp,μ,-Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s + 2*ω*Ep(p))/(2*p*q))^2))
		result = hquadrature(integrand1,p1,p2,reltol=1e-2,maxevals=10000)[1]

		integrand2(p) = -(s-4*m^2)*(numberF(temp,μ,-ω - Ep(p)) - numberF(temp,μ,-Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s + 2*ω*Ep(p))/(2*p*q))^2)) - (s-4*m^2)*(numberF(temp,μ,-ω + Ep(p)) - numberF(temp,μ,Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s - 2*ω*Ep(p))/(2*p*q))^2))
		result += hquadrature(integrand2,p2,param.Λ,reltol=1e-2,maxevals=10000)[1]
	end

	return result
end
function imagpart_σ_q(temp,μ,ω,q,m,param)
	β = 1/temp
	s = ω^2 - q^2
	Ep(p) = sqrt(p^2 + m^2)

	if s>=0
		if s<=4*m^2 || ω > 2*param.Λ
			result = 0.0
		else
			p1 = abs(-q + ω*sqrt(1 - 4*m^2/s))/2
			p2 = (q + ω*sqrt(1 - 4*m^2/s))/2
			integrand(p) = (s-4*m^2)*(numberF(temp,μ,-ω + Ep(p)) - numberF(temp,μ,Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s - 2*ω*Ep(p))/(2*p*q))^2))

			result = hquadrature(integrand,p1,p2,reltol=1e-2,maxevals=10000)[1]
		end
	else 
		p1 = abs(-q + ω*sqrt(1 - 4*m^2/s))/2
		p2 = (q + ω*sqrt(1 - 4*m^2/s))/2
		integrand1(p) = -(s-4*m^2)*(numberF(temp,μ,-ω - Ep(p)) - numberF(temp,μ,-Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s + 2*ω*Ep(p))/(2*p*q))^2))
		result = hquadrature(integrand1,p1,p2,reltol=1e-2,maxevals=10000)[1]

		integrand2(p) = -(s-4*m^2)*(numberF(temp,μ,-ω - Ep(p)) - numberF(temp,μ,-Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s + 2*ω*Ep(p))/(2*p*q))^2)) - (s-4*m^2)*(numberF(temp,μ,-ω + Ep(p)) - numberF(temp,μ,Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s - 2*ω*Ep(p))/(2*p*q))^2))
		result += hquadrature(integrand2,p2,param.Λ,reltol=1e-2,maxevals=10000)[1]
	end

	return result
end

"""
	realpart_σ_q(temp, μ, ω,param)

Returns the real part of the polarisation for σ₁ and σ₂ at external momentum q.
"""
function realpart_σ_q(temp, μ, ω, q, param)

	integrand(ν) = 2*ν*imagpart_σ_q(temp,μ,ν,q,param)*(PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2))/π
	int(ν) = integrand(1/(1 - ν))/(1 - ν)^2

	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]  + hquadrature(int,0,1,reltol=1e-2,maxevals=10000)[1]
end

function realpart_σ_q(temp,μ, ω, q, m, param)
	integrand(ν) = 2*ν*imagpart_σ_q(temp,μ,ν,q,m,param)*(PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2))
	int(ν) = integrand(1/(1 - ν))/(1 - ν)^2

	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]  + hquadrature(int,0,1,reltol=1e-2,maxevals=10000)[1]
end

@doc raw"""
	phasesc_σ_q(temp,μ,ω,param)

Returns the scattering part of the phase for σ₁ and σ₂ at external momentum q.

Also consider the equivalent function 

    phasesc_ϕ_q(temp,μ,ω,m,param)

where you supply the values of $\bar{\sigma}_1 = m$. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phasesc_σ_q(temp,μ, ω, q, param)
	repi = realpart_σ_q(temp,μ,ω,q,param)
	impi = imagpart_σ_q(temp,μ,ω,q,param)
	return angle(Complex(repi,-impi))
end

function phasesc_σ_q(temp,μ, ω, q, m, param)
	repi = realpart_σ_q(temp,μ,ω,q,m,param)
	impi = imagpart_σ_q(temp,μ,ω,q,m,param)

	return angle(Complex(repi,-impi))
end

@doc raw"""
	phaser_σ_q(temp,μ,ω,param)

Returns the resonant part of the phase for σ₁ and σ₂ at external momentum q.

Also consider the equivalent function 

    phaser_ϕ_q(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency 
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases
at different vaules of frequencies at a fixed temp and μ.
"""
function phaser_σ_q(temp,μ, ω, q, param)
	repi = realpart_σ_q(temp,μ,ω,q,param)
	impi = imagpart_σ_q(temp,μ,ω,q,param)
	Π00 = Π0_σ(temp,μ,param)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

function phaser_σ_q(temp,μ, ω, q, m, Π00, param)
	repi = realpart_σ_q(temp,μ,ω,q,m,param)
	impi = imagpart_σ_q(temp,μ,ω,q,m,param)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

@doc raw"""
	phase_σ_q(temp,μ,ω,param)

Returns all the phases for σ₁ and σ₂ at external momentum q in an array 

	[scattered phase, resonant phase, total phase]

Also consider the equivalent function 

    phase_σ_q(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phase_σ_q(temp,μ, ω, q, param)
	repi = realpart_σ_q(temp,μ,ω,q,param)
	impi = imagpart_σ_q(temp,μ,ω,q,param)
	Π00 = Π0_σ(temp,μ,param)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

function phase_σ_q(temp,μ, ω, q, m, Π00, param)
	repi = realpart_σ_q(temp,μ,ω,q,m,param)
	impi = imagpart_σ_q(temp,μ,ω,q,m,param)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end


export imagpart_ϕ_q, imagpart_σ_q, realpart_ϕ_q, realpart_σ_q, phasesc_ϕ_q, phasesc_σ_q, phaser_ϕ_q, phaser_σ_q, phase_ϕ_q ,phase_σ_q