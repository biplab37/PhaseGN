# This file contains code for solving the gap equation for the NJL model in 3D.

function energy(p, m, delta, mu)
    ep = sqrt(p^2 + m^2)
    return sqrt((ep + mu)^2 + delta^2)
end

## Normal Phase Δ_MF = 0
Enp(p, m, mu) = sqrt(p^2 + m^2) + mu
Enm(p, m, mu) = sqrt(p^2 + m^2) - mu
En(p, m) = sqrt(p^2 + m^2)

function integrand_m(T, μ, m, ω, param::Parameters3D)
    return p -> 4 * 6 * param.coupling.Gs * (p^2 / 2π^2) * m *
                (1 - numberF(T, μ, Enp(p, m, μ) + ω) - numberF(T, μ, Enm(p, m, μ) - ω)) / En(p, m)
end

function integrand_μ(T, μ, m, ω, param::Parameters3D)
    return p -> 4 * 6 * param.coupling.Gv * (p^2 / 2π^2) * m *
                (numberF(T, μ, En(p, m) - μ - ω) - numberF(T, μ, En(p, m) + μ + ω))
end

function gapeqns(T, μ, param::Parameters3D)
    function gap(m, ω)
        return [
            m - param.m0 - integrate(integrand_m(T, μ, m, ω, param), 0, param.cutoff),
            ω - integrate(integrand_μ(T, μ, m, ω, param), 0, param.cutoff),
        ]
    end
    return x -> gap(x[1], x[2])
end

"""
    massgap(T, μ, param::Parameters)
    massgap(trange::AbstractRange, μ, param::Parameters; initial_guess = [0.4, 0.0])
    massgap(T, μrange::AbstractRange, param::Parameters; initial_guess = [0.4, 0.0])
    massgap(trange::AbstractRange, μrange::AbstractRange, param::Parameters; initial_guess = [0.4, 0.0])

Solves the gap equation in normal phase for a given temperature and chemical potential.
Returns the quark mass and the ω condensate.
"""
function massgap(T, μ, param::Parameters3D)
    return nlsolve(gapeqns(T, μ, param), [0.4, 0.1]).zero
end

function massgap(trange::AbstractRange, μ, param::Parameters3D; initial_guess=[0.4, 0.0])
    result_m = zeros(length(trange))
    result_ω = zeros(length(trange))
    sol = initial_guess
    for (i, t) in enumerate(trange)
        sol = nlsolve(gapeqns(t, μ, param), sol).zero
        result_m[i] = sol[1]
        result_ω[i] = sol[2]
    end
    return result_m, result_ω
end

function massgap(T, μrange::AbstractRange, param::Parameters3D; initial_guess=[0.4, 0.0])
    result_m = zeros(length(μrange))
    result_ω = zeros(length(μrange))
    sol = initial_guess
    for (i, μ) in enumerate(μrange)
        sol = nlsolve(gapeqns(T, μ, param), sol).zero
        result_m[i] = sol[1]
        result_ω[i] = sol[2]
    end
    return result_m, result_ω
end

function massgap(trange::AbstractRange, μrange::AbstractRange, param::Parameters3D; initial_guess=[0.4, 0.0])
    if length(trange) != length(μrange)
        throw(DimensionMismatch("temperature list and the chemical potential list must have the same length"))
    end
    result_m = zeros(length(trange))
    result_ω = zeros(length(trange))
    sol = initial_guess
    for (i, (t, μ)) in enumerate(zip(trange, μrange))
        sol = nlsolve(gapeqns(t, μ, param), sol).zero
        result_m[i] = sol[1]
        result_ω[i] = sol[2]
    end
    return result_m, result_ω
end

export massgap
