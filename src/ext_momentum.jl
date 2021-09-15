# This file contains the code for the case with non zero external momentum

#-------------------------------------------------------------#
#|                      Pseudo Scalar(ϕ)                     |#
#-------------------------------------------------------------#


"""
	imagpart_ϕ_q(temp,μ,ω,q,κ)

Returns the imaginary part of the polarisation function at finite external momentum q for the pseudo scalar channel.
"""
function imagpart_ϕ_q(temp,μ,ω,q,κ)
	β = 1/temp
	m = σ1(temp,μ,κ)
	s = ω^2 - q^2
	Ep(p) = sqrt(p^2 + m^2)

	if s>=0
		if s<=4*m^2 || ω > 2*Λ[1]
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

		result += hquadrature(integrand2,p2,Λ[1],reltol=1e-2,maxevals=10000)[1]
	end

	return result
end
function imagpart_ϕ_q(temp,μ,ω,q,κ,m)
	β = 1/temp
	s = ω^2 - q^2
	Ep(p) = sqrt(p^2 + m^2)

	if s>=0
		if s<=4*m^2 || ω > 2*Λ[1]
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
	realpart_ϕ_q(temp, μ, ω, q, κ)

Returns the real part of the polarisation function at finite external momentum for the pseudo scalar channel.
Uses `imagpart_ϕ_q` and Kramer-Cronig Relation for the calculation
"""
function realpart_ϕ_q(temp, μ, ω, q, κ)
	
	integrand(ν) = 2*ν*imagpart_ϕ_q(temp,μ,ν,q,κ)*(PrincipleValue(ν^2 - ω^2) - PrincipleValue(ν^2))/π
	int1(ν) = integrand(1/(1 - ν))/(1 - ν)^2

	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]  + hquadrature(int1,0,1,reltol=1e-2,maxevals=10000)[1]
end

function realpart_ϕ_q(temp, μ, ω, q, κ, m)
	
	integrand(ν) = 2*ν*imagpart_ϕ_q(temp,μ,ν,q,κ,m)*(PrincipleValue(ν^2 - ω^2) - PrincipleValue(ν^2))/π
	int1(ν) = integrand(1/(1 - ν))/(1 - ν)^2

	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]  + hquadrature(int1,0,1,reltol=1e-2,maxevals=10000)[1]
end


function phasesc_ϕ_q(temp,μ,ω,q,κ)
	repi = realpart_ϕ_q(temp,μ,ω,q,κ)
	impi = imagpart_ϕ_q(temp,μ,ω,q,κ)
	return angle(Complex(repi,-impi))
end
function phasesc_ϕ_q(temp,μ,ω,q,κ,m)
	repi = realpart_ϕ_q(temp,μ,ω,q,κ,m)
	impi = imagpart_ϕ_q(temp,μ,ω,q,κ,m)
	
	return angle(Complex(repi,-impi))
end

function phaser_ϕ_q(temp,μ,ω,q,κ)
    repi = realpart_ϕ_q(temp,μ,ω,q,κ)
	impi = imagpart_ϕ_q(temp,μ,ω,q,κ)
	Π00 = Π0_ϕ(temp,μ,κ)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end
function phaser_ϕ_q(temp,μ,ω,q,κ,m,Π00)
    repi = realpart_ϕ_q(temp,μ,ω,q,κ,m)
	impi = imagpart_ϕ_q(temp,μ,ω,q,κ,m)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

function phase_ϕ_q(temp,μ,ω,q,κ)
    repi = realpart_ϕ_q(temp,μ,ω,q,κ)
	impi = imagpart_ϕ_q(temp,μ,ω,q,κ)
	Π00 = Π0_ϕ(temp,μ,κ)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

function phase_ϕ_q(temp,μ,ω,q,κ,m,Π00)
    repi = realpart_ϕ_q(temp,μ,ω,q,κ,m)
	impi = imagpart_ϕ_q(temp,μ,ω,q,κ,m)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

#-------------------------------------------------------------#
#|                         Scalar(σ)                         |#
#-------------------------------------------------------------#

function imagpart_σ_q(temp,μ,ω,q,κ)
	β = 1/temp
	m = σ1(temp,μ,κ)
	s = ω^2 - q^2
	Ep(p) = sqrt(p^2 + m^2)

	if s>=0
		if s<=4*m^2 || ω > 2*Λ[1]
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
		result += hquadrature(integrand2,p2,Λ[1],reltol=1e-2,maxevals=10000)[1]
	end

	return result
end
function imagpart_σ_q(temp,μ,ω,q,κ,m)
	β = 1/temp
	s = ω^2 - q^2
	Ep(p) = sqrt(p^2 + m^2)

	if s>=0
		if s<=4*m^2 || ω > 2*Λ[1]
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
		result += hquadrature(integrand2,p2,Λ[1],reltol=1e-2,maxevals=10000)[1]
	end

	return result
end

function realpart_σ_q(temp, μ, ω, q, κ)

	integrand(ν) = 2*ν*imagpart_σ_q(temp,μ,ν,q,κ)*(PrincipleValue(ν^2 - ω^2) - PrincipleValue(ν^2))/π
	int(ν) = integrand(1/(1 - ν))/(1 - ν)^2

	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]  + hquadrature(int,0,1,reltol=1e-2,maxevals=10000)[1]
end

function realpart_σ_q(temp,μ, ω, q, κ, m)
	integrand(ν) = 2*ν*imagpart_σ_q(temp,μ,ν,q,κ,m)*(PrincipleValue(ν^2 - ω^2) - PrincipleValue(ν^2))
	int(ν) = integrand(1/(1 - ν))/(1 - ν)^2

	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]  + hquadrature(int,0,1,reltol=1e-2,maxevals=10000)[1]
end

function phasesc_σ_q(temp,μ, ω, q, κ)
	repi = realpart_σ_q(temp,μ,ω,q,κ)
	impi = imagpart_σ_q(temp,μ,ω,q,κ)
	return angle(Complex(repi,-impi))
end

function phasesc_σ_q(temp,μ, ω, q, κ, m)
	repi = realpart_σ_q(temp,μ,ω,q,κ,m)
	impi = imagpart_σ_q(temp,μ,ω,q,κ,m)

	return angle(Complex(repi,-impi))
end

function phaser_σ_q(temp,μ, ω, q, κ)
	repi = realpart_σ_q(temp,μ,ω,q,κ)
	impi = imagpart_σ_q(temp,μ,ω,q,κ)
	Π00 = Π0_σ(temp,μ,κ)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

function phaser_σ_q(temp,μ, ω, q, κ, m, Π00)
	repi = realpart_σ_q(temp,μ,ω,q,κ,m)
	impi = imagpart_σ_q(temp,μ,ω,q,κ,m)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

function phase_σ_q(temp,μ, ω, q, κ)
	repi = realpart_σ_q(temp,μ,ω,q,κ)
	impi = imagpart_σ_q(temp,μ,ω,q,κ)
	Π00 = Π0_σ(temp,μ,κ)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

function phase_σ_q(temp,μ, ω, q, κ, m, Π00)
	repi = realpart_σ_q(temp,μ,ω,q,κ,m)
	impi = imagpart_σ_q(temp,μ,ω,q,κ,m)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end


export imagpart_ϕ_q, imagpart_σ_q, realpart_ϕ_q, realpart_σ_q, phasesc_ϕ_q, phasesc_σ_q, phaser_ϕ_q, phaser_σ_q, phase_ϕ_q ,phase_σ_q