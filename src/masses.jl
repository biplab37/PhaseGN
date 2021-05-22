@doc raw"""
    σ1(temp,μ, κ=0.00, M=1.0)

Returns the value of $\bar{\sigma}_1$ by solving the gap equation.
"""
function σ1(temp,μ, κ=0.00, M=1.0)
    β = 1/temp
	
	## The gap equation
    σ₁(σ) = (1 + 2*cosh(β*μ)*exp(-β*σ) + exp(-2*β*σ) - exp(β*(M - σ + (π*κ/σ))))
	
    result = fzero(σ₁,1e-4,4)
	
    ## For the case when there is no zero in the interval
	if (4.0-result<1e-3 || result == 1e-4)
        result = 0.0
    end
	
    return result
end

function gE2integ(temp, μ, κ, p,ME)
    β = 1/temp
    m = σ1(temp,μ,κ)
    Ep = sqrt(m^2 + p^2)

    return p*(1 - 1/(1+exp(β*(Ep+μ))) - 1/(1+exp(β*(Ep-μ))))*PrincipleValue(2*Ep - ME)/(Ep*(2*Ep + ME)*pi)
end

"""
	gE2(temp,μ,κ,ME)

Returns the square of the induced Fermion-exciton coupling constant.
"""
function gE2(temp,μ, κ,ME)
    func(p) = gE2integ(temp,μ, κ, p,ME)
    return 1/hquadrature(func,0,Λ[1],maxevals=10000)[1]
end

"""
	M_σ(temp,μ,κ)

Returns the exciton mass of scalar channels σ₁ and σ₂
"""
function M_σ(temp,μ,κ)
    m = σ1(temp,μ,κ)
    func(ME) = ME^2 - 4*m^2 - κ*gE2(temp,μ,κ,ME)/m
    result = fzero(func,0.0,25)
    if (result == 20 || result == 0.0)
        print("mass not found")
    else
        return result
    end
end

"""
	M_ϕ(temp,μ,κ)

Returns the exciton mass of scalar channels ϕ₁ and ϕ₂
"""
function M_ϕ(temp,μ,κ)
    m = σ1(temp,μ,κ)
    func(ME) = ME^2 - κ*gE2(temp,μ,κ,ME)/m
    result = fzero(func,0.0,25)
    if (result == 20 || result == 0.0)
        print("mass not found")
    else
        return result
    end
end

export σ1, M_ϕ, M_σ
