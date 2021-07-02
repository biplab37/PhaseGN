# This file contains the code for the case with non zero external momentum

#-------------------------------------------------------------#
#|                      Pseudo Scalar(ϕ)                     |#
#-------------------------------------------------------------#


"""
	imagpart_ϕ(temp,μ,ω,q,κ)

Returns the imaginary part of the polarisation function at finite external momentum q for the pseudo scalar channel.
"""
function imagpart_ϕ(temp,μ,ω,q,κ)
	β = 1/temp
	m = σ1(temp,μ,κ)
	s = ω^2 - q^2
	Ep(p) = sqrt(p^2 + m^2)

	if s>=0
		if s<=4*m^2 || ω > 2*Λ
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

		integrand2(p) = -s*(numberF(temp,μ,-ω - Ep(p)) - numberF(temp,μ,-Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s + 2*ω*Ep(p))/(2*p*q))^2)) + s*(numberF(temp,μ,-ω + Ep(p)) - numberF(temp,μ,Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s - 2*ω*Ep(p))/(2*p*q))^2))
		int1(p) = integrand2(1/(1 - p))/(1 - p)^2

		result += hquadrature(int1,1 - 1/p2,1,reltol=1e-2,maxevals=10000)[1]
	end

	return result
end
function imagpart_ϕ(temp,μ,ω,q,κ,m)
	β = 1/temp
	s = ω^2 - q^2
	Ep(p) = sqrt(p^2 + m^2)

	if s>=0
		if s<=4*m^2 || ω > 2*Λ
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

		integrand2(p) = -s*(numberF(temp,μ,-ω - Ep(p)) - numberF(temp,μ,-Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s + 2*ω*Ep(p))/(2*p*q))^2)) + s*(numberF(temp,μ,-ω + Ep(p)) - numberF(temp,μ,Ep(p)))/(8*π*q*Ep(p)*sqrt(1 - ((s - 2*ω*Ep(p))/(2*p*q))^2))
		int1(p) = integrand2(1/(1 - p))/(1 - p)^2

		result += hquadrature(int1,1 - 1/p2,1,reltol=1e-2,maxevals=10000)[1]
	end

	return result
end


"""
	realpart_ϕ(temp, μ, ω, q, κ)

Returns the real part of the polarisation function at finite external momentum for the pseudo scalar channel.
Uses `imagpart_ϕ` and Kramer-Cronig Relation for the calculation
"""
function realpart_ϕ(temp, μ, ω, q, κ)
	
	integrand(ν) = 2*ν*imagpart_ϕ(temp,μ,ν,q,κ)*(PrincipleValue(ν^2 - ω^2) - PrincipleValue(ν^2))/π
	int1(ν) = integrand(1/(1 - ν))/(1 - ν)^2

	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]  + hquadrature(int1,0,1,reltol=1e-2,maxevals=10000)[1]
end

function realpart_ϕ(temp, μ, ω, q, κ, m)
	
	integrand(ν) = 2*ν*imagpart_ϕ(temp,μ,ν,q,κ,m)*(PrincipleValue(ν^2 - ω^2) - PrincipleValue(ν^2))/π
	int1(ν) = integrand(1/(1 - ν))/(1 - ν)^2

	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]  + hquadrature(int1,0,1,reltol=1e-2,maxevals=10000)[1]
end


function phasesc_ϕ(temp,μ,ω,q,κ)
	repi = realpart_ϕ(temp,μ,ω,q,κ)
	impi = imagpart_ϕ(temp,μ,ω,q,κ)
	return angle(Complex(repi,-impi))
end
function phasesc_ϕ(temp,μ,ω,q,κ,m)
	repi = realpart_ϕ(temp,μ,ω,q,κ,m)
	impi = imagpart_ϕ(temp,μ,ω,q,κ,m)
	
	return angle(Complex(repi,-impi))
end

function phaser_ϕ(temp,μ,ω,q,κ)
    repi = realpart_ϕ(temp,μ,ω,q,κ)
	impi = imagpart_ϕ(temp,μ,ω,q,κ)
	Π00 = Π0_ϕ(temp,μ,κ)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end
function phaser_ϕ(temp,μ,ω,q,κ,m,Π00)
    repi = realpart_ϕ(temp,μ,ω,q,κ,m)
	impi = imagpart_ϕ(temp,μ,ω,q,κ,m)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

function phase_ϕ(temp,μ,ω,q,κ)
    repi = realpart_ϕ(temp,μ,ω,q,κ)
	impi = imagpart_ϕ(temp,μ,ω,q,κ)
	Π00 = Π0_ϕ(temp,μ,κ)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

function phase_ϕ(temp,μ,ω,q,κ,m,Π00)
    repi = realpart_ϕ(temp,μ,ω,q,κ,m)
	impi = imagpart_ϕ(temp,μ,ω,q,κ,m)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

#-------------------------------------------------------------#
#|                         Scalar(σ)                         |#
#-------------------------------------------------------------#

