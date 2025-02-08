# Correct regularization
using PhaseGN
using Plots

function omega(s, temp, μ, param)
    term_ind = abs(s)^3 / (3 * π) - s^2 * param.M / (2 * π) #- param.κ * s - s^4 / (8 * π * param.Λ)
    σ₁ = abs(s)
    term_dep = temp^3 * (σ₁ * (PhaseGN.reli2(-exp(-(σ₁ - μ) / temp)) + PhaseGN.reli2(-exp(-(σ₁ + μ) / temp))) / temp + PhaseGN.reli3(-exp(-(σ₁ - μ) / temp)) + PhaseGN.reli3(-exp(-(σ₁ + μ) / temp))) / π
    correction_term = temp^3 * (-PhaseGN.reli3(-exp(μ / temp)) - PhaseGN.reli3(-exp(-μ / temp))) / π
    return term_dep + term_ind
end

srange = -1.5:0.01:1.5

param = Parameters(Λ=5.0, κ=0.046)

ome = omega.(srange, 0.01, 1.01, param)

plot(srange, ome)

function second_derivative(f, x; dx=1e-2)
    return (f(x + dx) + f(x - dx) - 2 * f(x)) / dx^2
end
function sg_filter(f, x; dx=5e-2)
    return (-1 * f(x - 2 * dx) + 16 * f(x - dx) - 30 * f(x) + 16 * f(x + dx) - 1 * f(x + 2 * dx)) / (12 * dx^2)
end
function first_derivative(f, x; dx=1e-2)
    return (-2 * f(x - 2 * dx) - f(x + dx) + f(x + dx) + 2 * f(x + 2 * dx)) / (10 * dx)
end

f1(x) = 1 / 10 * ((1 + x)^2 / x)
xrange = 0.5:0.01:20.0
plot(xrange, [f1.(xrange) log.(f1.(20 .* xrange))])

### Plotting the Tricritical point

function criticalLine(T, param=Parameters())
    return (T < (param.M / (2 * log(2)))) ? T * acosh(0.5 * exp(1 / T) - 1.0) : 0.0
end

function second_line(T, param=Parameters())
    return 2 * T * atanh(sqrt(1 - 2 * T / param.Λ))
end

mus1 = criticalLine.(0.01:0.01:0.8, param)
mus2 = second_line.(0.01:0.01:0.8, Parameters(Λ=2.0))

plot(0.01:0.01:0.8, [mus1 mus2])