function imagpart_ϕ(temp, μ, ω,κ)
	m = σ1(temp,μ,κ)
	Nσ = ω^2 - 4*m^2
	if Nσ > 0 && abs(ω) < 2*sqrt(Λ[1]^2 + m^2)
		return ω*(1 - numberF(temp,-μ,ω/2) - numberF(temp,μ,ω/2))/(4)
	else
		return 0.0
	end
end

function realpart_ϕ(temp,μ,ω,κ)
	m = σ1(temp,μ,κ)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = - p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipleValue(ω - 2*Ep(p)) - PrincipleValue(ω + 2*Ep(p)) + PrincipleValue(Ep(p)))/π
    return hquadrature(integrand,0.0,Λ[1],reltol=1e-3,maxevals=10000)[1]
end

function fullrealpart_ϕ(temp,μ,ω,κ)
	m = σ1(temp,μ,κ)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = 1/π + p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipleValue(ω - 2*Ep(p)) - PrincipleValue(ω + 2*Ep(p)))/π
    return  hquadrature(integrand,0.0,Λ[1],reltol=1e-3,maxevals=10000)[1]
end

function Π0_ϕ(temp,μ,κ)
	m = σ1(temp,μ,κ)
	Ep(p) = sqrt(p^2+m^2)
	integrand(p) = 1/π - p*(1)*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))/(π*Ep(p))
	return hquadrature(integrand,0.0,Λ[1],reltol=1e-3,maxevals=10000)[1]
end

function phasesc_ϕ(temp,μ,ω,κ)
	repi = realpart_ϕ(temp,μ,ω,κ)
	impi = imagpart_ϕ(temp,μ,ω,κ)
	return angle(Complex(repi,-impi))
end

function phaser_ϕ(temp,μ,ω,κ)
    repi = realpart_ϕ(temp,μ,ω,κ)
	impi = imagpart_ϕ(temp,μ,ω,κ)
	Π00 = Π0_ϕ(temp,μ,κ)

	return -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))
end

function phase_ϕ(temp,μ,ω,κ)
    repi = realpart_ϕ(temp,μ,ω,κ)
	impi = imagpart_ϕ(temp,μ,ω,κ)
	Π00 = Π0_ϕ(temp,μ,κ)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

function phase_ϕ(temp,μ,ω,κ,Π00)
    repi = realpart_ϕ(temp,μ,ω,κ)
	impi = imagpart_ϕ(temp,μ,ω,κ)

	phasesc = angle(Complex(repi,-impi))
	phaser = -angle(Complex(repi - (repi^2 + impi^2)/Π00,-impi))

	return [phasesc, phaser, phasesc + phaser]
end

export imagpart_ϕ, realpart_ϕ, fullrealpart_ϕ, phasesc_ϕ, phaser_ϕ, phase_ϕ