@doc raw"""
    σ1(temp,μ, p::Parameters)

Returns the value of $\bar{\sigma}_1$ by solving the gap equation.
"""
function σ1(temp::Float64, μ, param::Parameters; guess=0.5)
    β = 1 / temp

    ## The gap equation
    σ₁(σ) = (1 + 2 * cosh(β * μ) * exp(-β * σ) + exp(-2 * β * σ) - exp(β * (param.M - σ + (π * param.κ / σ))))

    result = bisection(σ₁, 1e-4, 4)

    ## For the case when there is no zero in the interval
    if (4.0 - result < 1e-3 || result == 1e-4)
        result = 0.0
    end

    return result
end

function σ1(temp, μ; guess=0.5)
    param = Parameters()
    @warn("No parameters given, using default parameters: κ = $(param.κ), M = $(param.M)")
    σ1(temp, μ, param, guess=guess)
end

function σ1(trange::AbstractRange, μ, param::Parameters; guess=0.5)
    σlist = zeros(length(trange))

    ## The gap equation
    σ₁(β, σ) = (1 + 2 * cosh(β * μ) * exp(-β * σ) + exp(-2 * β * σ) - exp(β * (param.M - σ + (π * param.κ / σ))))

    for (i, t) in enumerate(trange)
        guess = fzero(σ -> σ₁(1 / t, σ), guess)
        σlist[i] = guess
    end

    return σlist
end

function critical_line(t::Number, param::Parameters)
    determinant = (exp(param.M / t) - 2)^2 - 4.0
    if determinant > 0.0
        return t * log(((exp(param.M / t) - 2) + sqrt(determinant)) / 2)
    else
        return 0.0
    end
end


function gE2integ(temp, μ, p, ME, param::Parameters)
    β = 1 / temp
    m = σ1(temp, μ, param)
    Ep = sqrt(m^2 + p^2)

    return p * (1 - 1 / (1 + exp(β * (Ep + μ))) - 1 / (1 + exp(β * (Ep - μ)))) * PrincipalValue(2 * Ep - ME) / (Ep * (2 * Ep + ME) * π)
end

"""
	gE2(temp,μ,ME,param::Parameters)

Returns the square of the induced Fermion-exciton coupling constant.
"""
function gE2(temp, μ, ME, param::Parameters)
    func(p) = gE2integ(temp, μ, p, ME, param)
    return 1 / integrate(func, 0, param.Λ)
end

"""
	M_σ(temp,μ,param::Parameters)

Returns the exciton mass of scalar channels σ₁ and σ₂
"""
function M_sigma(temp, μ, param::Parameters)
    m = σ1(temp, μ, param)
    func(ME) = ME^2 - 4 * m^2 - param.κ * gE2(temp, μ, ME, param) / m
    result = bisection(func, 0.0, 5)
    if (result ≈ 5 || result == 0.0)
        print("mass not found")
    else
        return result
    end
end

"""
	M_ϕ(temp,μ,param::Parameters)

Returns the exciton mass of scalar channels ϕ₁ and ϕ₂
"""
function M_phi(temp, μ, param::Parameters)
    m = σ1(temp, μ, param)
    func(ME) = ME^2 - param.κ * gE2(temp, μ, ME, param) / m
    result = bisection(func, 0.0, 5)
    if (result == 5 || result == 0.0)
        print("mass not found")
    else
        return result
    end
end

export σ1, M_phi, M_sigma
