function realpartI2(temp,μ,κ)
    m = σ1(temp,μ,κ)
    M = M_ϕ(temp,μ,κ)
    Ep(p) = sqrt(p^2+m^2)
    integrand(p) = p*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))*(PrincipleValue(2*Ep(p) - M) + PrincipleValue(M + 2*Ep(p)))/(8*π*Ep(p)^2)
    return hquadrature(integrand,0.0,Λ[1],reltol=1e-2,maxevals=10000)[1]
end

function imagpartI2(temp, μ,κ)
    m = σ1(temp,μ,κ)
    M = M_ϕ(temp,μ,κ)
    Nσ = M^2 - 4*m^2
    if Nσ > 0 && abs(M) < 2*sqrt(Λ[1]^2 + m^2)
        return (1 - numberF(temp,-μ,M/2) - numberF(temp,μ,M/2))/(8*M)
    else
        return 0.0
    end
end

function I1(temp,μ,κ)
    m = σ1(temp,μ,κ)
    Ep(p) = sqrt(p^2+m^2)
    integrand(p) = 1/(2*π) - 2*p*(1 - numberF(temp,-μ,Ep(p)) - numberF(temp,μ,Ep(p)))/(4*π*Ep(p))
    return -1/(2*π) + hquadrature(integrand,0.0,Λ[1],reltol=1e-2,maxevals=10000)[1]
end

function Γϕ(temp,μ,κ)
    M = M_ϕ(temp,μ,κ)
    imp = imagpartI2(temp,μ,κ)
    rep = realpartI2(temp,μ,κ)
    i1 = I1(temp,μ,κ)

    Γ = (i1*imp/((imp^2+rep^2)*M))
    return Γ
end

function Γϕ2(temp,μ,κ)
	m = σ1(temp,μ,κ)
    M = M_ϕ(temp,μ,κ)
    imp = imagpartI2(temp,μ,κ)
    rep = realpartI2(temp,μ,κ)
    i1 = I1(temp,μ,κ)
    x = M^2 - (i1*rep/((imp^2+rep^2)))

    if x>=0
    	return 2*sqrt(x)
    else
    	return 0.0
    end
end

export Γϕ, Γϕ2