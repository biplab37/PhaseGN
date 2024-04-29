"""
	pressure(func,temp,μ,param::Parameters)

Returns the pressure corrosponding to the phase function `func(temp,μ,ω,param)` 
"""
function pressure(func::Function, temp, μ, param::Parameters)
    integrand(s) = -(log(exp(sqrt(s) / temp) - 1) - sqrt(s) / temp) / (2 * π^2 * temp^2) * func(temp, μ, sqrt(s), param)
    int1(s) = integrand(1 / (1 - s)) / (1 - s)^2
    return integrate(integrand, 0, 1) + integrate(int1, 0, 1)
end

function pressure(func::Function, temp, μ, m, param::Parameters)
    integrand(s) = -(log(exp(sqrt(s) / temp) - 1) - sqrt(s) / temp) / (2 * π^2 * temp^2) * func(temp, μ, sqrt(s), m, param)
    int1(s) = integrand(1 / (1 - s)) / (1 - s)^2
    return integrate(integrand, 0, 1)
end

function pressure(func::Function, temp, μ, m, M, Γ, param::Parameters)
    integrand(s) = -(log(exp(sqrt(s) / temp) - 1) - sqrt(s) / temp) / (2 * π^2 * temp^2) * func(temp, μ, sqrt(s)m, M, Γ, param)
    int1(s) = integrand(1 / (1 - s)) / (1 - s)^2
    return integrate(integrand, 0, 1)
end

function pressure_fl_q(func::Function, temp, μ, q, param::Parameters)
    integrand(s) = (1 + 2 / (exp(sqrt(s + q^2) / temp) - 1)) * func(sqrt(s + q^2), temp, μ, q, param) / (2π)
    return integrate(integrand, -q^2, param.Λ^2)
end

export pressure
