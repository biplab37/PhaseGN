@doc raw"""
    σ1(temp,μ, p::Parameters)

Returns the value of $\bar{\sigma}_1$ by solving the gap equation.
"""
function σ1(temp,μ, param::Parameters)
    β = 1/temp
	
	## The gap equation
    σ₁(σ) = (1 + 2*cosh(β*μ)*exp(-β*σ) + exp(-2*β*σ) - exp(β*(param.M - σ + (π*param.κ/σ))))
	
    result = bisection(σ₁,1e-4,4)
	
    ## For the case when there is no zero in the interval
	if (4.0-result<1e-3 || result == 1e-4)
        result = 0.0
    end
	
    return result
end

function σ1(temp,μ)
    param = Parameters()
    @warn("No parameters given, using default parameters: κ = 0.0, M = 1.0")
    β = 1/temp
	
	## The gap equation
    σ₁(σ) = (1 + 2*cosh(β*μ)*exp(-β*σ) + exp(-2*β*σ) - exp(β*(param.M - σ + (π*param.κ/σ))))
	
    result = bisection(σ₁,1e-4,4)
	
    ## For the case when there is no zero in the interval
	if (4.0-result<1e-3 || result == 1e-4)
        result = 0.0
    end
	
    return result
end


function gE2integ(temp, μ, p,ME,param::Parameters)
    β = 1/temp
    m = σ1(temp,μ,param)
    Ep = sqrt(m^2 + p^2)

    return p*(1 - 1/(1+exp(β*(Ep+μ))) - 1/(1+exp(β*(Ep-μ))))*PrincipleValue(2*Ep - ME)/(Ep*(2*Ep + ME)*π)
end

"""
	gE2(temp,μ,ME,param::Parameters)

Returns the square of the induced Fermion-exciton coupling constant.
"""
function gE2(temp,μ,ME,param::Parameters)
    func(p) = gE2integ(temp,μ,p,ME, param)
    return 1/hquadrature(func,0,param.Λ,reltol=1e-3,maxevals=10000)[1]
end

"""
	M_σ(temp,μ,κ)

Returns the exciton mass of scalar channels σ₁ and σ₂
"""
function M_σ(temp,μ,param::Parameters)
    m = σ1(temp,μ,param)
    func(ME) = ME^2 - 4*m^2 - κ*gE2(temp,μ,ME,param)/m
    result = fzero(func,0.0,25)
    if (result ≈ 25 || result == 0.0)
        print("mass not found")
    else
        return result
    end
end

"""
	M_ϕ(temp,μ,κ)

Returns the exciton mass of scalar channels ϕ₁ and ϕ₂
"""
function M_ϕ(temp,μ,param::Parameters)
    m = σ1(temp,μ,param)
    func(ME) = ME^2 - κ*gE2(temp,μ,ME,param)/m
    result = fzero(func,0.0,25)
    if (result == 25 || result == 0.0)
        print("mass not found")
    else
        return result
    end
end

export σ1, M_ϕ, M_σ
