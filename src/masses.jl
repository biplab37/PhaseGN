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
    result = bisection(func, 0.0, param.Λ)
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
    result = bisection(func, 0.0, param.Λ)
    if (result == param.Λ || result == 0.0)
        print("mass not found")
    else
        return result
    end
end

function M_phi_bounded(temp, mu, param::Parameters)
    m = mass_k(temp, mu, 0.0, param)
    func(ME) = Π0_phi(temp, mu, param) - realpart(imagpart_phi_q_refactored_m, ME, temp, mu, 0.0, m, param)
    return bisection(func, 0.0, param.M)
end

function mass_k(temp, μ, k, param)
    β = 1 / temp
    integrand(m) = m * (m - param.M) - π * param.κ + angular_integral(temp, μ, k, m, param) / 2π
    return bisection(integrand, 0.0, 1.5 * param.M)
end

function angular_integral(temp, μ, k, m, param)
    β = 1 / temp
    f(sign, x) = numberF(temp, sign * μ, x)
    ksq(p) = p^2 + k^2 / 4 + m^2
    integrand(p, θ) = p * m * (f(+1.0, sqrt(ksq(p) + p * k * cos(θ))) + f(-1.0, sqrt(ksq(p) - p * k * cos(θ)))) / sqrt(p^2 + m^2)
    return integrate(x -> integrand(x[1], x[2]), [0.0, 0.0], [param.Λ, 2π])
end

"""
    find_kappa(cutoff; mp=0.1)

This function finds the value of kappa such that the pseudo scalar mass equals to mp.
"""
function find_kappa(cutoff; mp=0.1)
    func(kappa) = M_phi_bounded(0.01, 0.0, Parameters(Λ=cutoff, κ=kappa)) - mp
    sol = bisection(func, 0.0, 0.1)
    if abs(sol) < 1e-4
        return bisection(func, 0.0, 0.5)
    else
        return sol
    end
end

export σ1, M_phi, M_sigma, mass_k
