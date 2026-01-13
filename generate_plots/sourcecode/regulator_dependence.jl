using PhaseGN, Plots, LaTeXStrings, Roots
using PGFPlotsX
function second_term(T, mu)
    return -1 / (2 * π) + T * (log(1 + exp(-mu / T)) + log(1 + exp(mu / T))) / (2 * π)
end

function fourth_term(T, mu, Λ)
    return sech(mu / (2 * T))^2 / (16 * π * T) - 1 / (8 * π * Λ)
end

function sixth_term(T, mu)
    return (2 * sech(mu / (2 * T))^2 - 3 * sech(mu / (2 * T))^4) / (576 * π * T^3)
end

function crit_line(T)
    return PhaseGN.bisection(x -> second_term(T, x), 0.0, 1.5)
end

function four_line(T, cutoff)
    return PhaseGN.bisection(x -> fourth_term(T, x, cutoff), 0.0, 2.0)
end

trange = 0.01:0.01:0.8
plot(crit_line.(trange), trange)
plot!(four_line.(trange, 5.0), trange)

function zero_line(T, mu, lambda)
    a2 = second_term(T, mu)
    a4 = fourth_term(T, mu, lambda)
    a6 = sixth_term(T, m)

end

function tcp(Λ)
    return PhaseGN.bisection(x -> fourth_term(x, crit_line(x), Λ), 0.01, 0.8)
end

lambda_range = [2, 3, 5, 10, 100]
    tcps = [tcp(l) for l in lambda_range]

writedlm("fig_12_b_TCpoints.csv", hcat(lambda_range, tcps), ',')

scatter!(crit_line.(tcps), tcps)

lambdarange = 2 .^ (1:0.01:6)
mvalues = [x * (1 - sqrt(1 - 2 / x)) for x in lambdarange]
p1 = @pgf Axis(
    {
        xmode = "log",
        ylabel = L"m(T=0, \mu=0)/M",
        xlabel = L"\Lambda/M",
        xmin = 1.9,
        xmax = 2^6,
    },
    Plot(
        {
            no_marks,
        },
        Table(lambdarange, mvalues)
    ),
    HLine(
        {
            style = "dashed",
            opacity = 0.3,
        },
        1.0
    )
)

pgfsave("mass_lambda_dependence.pdf", p1)

## phase structure

function discriminant_cubic(T, mu, kappa)
    return
end

function cardano(a, b, kappa)
    # ax^3 + bx - kappa = 0
    if a == 0
        if b == 0
            @warn "No solution!"
        end
        return kappa / b
    end
    if b == 0
        return (kappa / a)^(1 / 3)
    end
    p = b / a
    q = -kappa / a

    Δ = q^2 / 4 + p^3 / 27
    if Δ < 0
        Δ = Complex(Δ)
    end

    return Complex(-q / 2 + sqrt(Δ))^(1 / 3) + Complex(-q / 2 - sqrt(Δ))^(1 / 3)
end
using Polynomials

# function find_roots(T, mu, kappa, lambda)
#     a1 = -kappa
#     a2 = second_term(T, mu)
#     a4 = fourth_term(T, mu, lambda)
#     a6 = sixth_term(T, mu)
#     if a4<0
#         pol = Polynomial([0.0, a1, a2, 0.0, a4, 0.0, a6])
#     else
#         pol = Polynomial([0.0, a1, a2, 0.0, a4])
#     end
#     derpol = derivative(pol)
#     rts = roots(derpol)
#     vals = pol.(rts)
#     return rts, vals
# end
function find_roots(T, mu, kappa, lambda)
    a1 = -kappa
    a2 = second_term(T, mu)
    a4 = fourth_term(T, mu, lambda)
    a6 = sixth_term(T, mu)
    if a4 < 0
        pol = Polynomial([0.0, a1, a2, 0.0, a4, 0.0, a6])
    else
        pol = Polynomial([0.0, a1, a2, 0.0, a4])
    end
    derpol = derivative(pol)
    rts = find_zeros(x -> derpol(x), (0.0, 2.1))
    vals = pol.(rts)
    return rts, vals
end

find_roots(0.14, 0.99, 0.01, 10.0)


function discriminant(T, mu, kappa, lambda)
    a1 = -kappa
    a2 = second_term(T, mu)
    a4 = fourth_term(T, mu, lambda)
    return -(4 * (a2 / (2 * a4))^3 + 27 * (a1 / (4 * a4))^3)
end

discriminant(0.8, 0.0, 0.046, 5.0)

trange = 0.01:0.01:0.8
murange = 0.0:0.01:1.5

plotly()
heatmap(murange, trange, (x, y) -> sign(fourth_term(y, x, 5.0)))

function omega(m, temp, μ, kappa, lambda, norm=false)
    σ₁ = abs(m)
    temp_independent = -(m^2 / 2π) + (abs(m)^3 / 3π) - kappa * m - m^4 / (8 * π * lambda)
    temp_dependent = temp^3 * (σ₁ * (PhaseGN.reli2(-exp(-(σ₁ - μ) / temp)) + PhaseGN.reli2(-exp(-(σ₁ + μ) / temp))) / temp + PhaseGN.reli3(-exp(-(σ₁ - μ) / temp)) + PhaseGN.reli3(-exp(-(σ₁ + μ) / temp))) / π
    vac_term = temp^3 * (PhaseGN.reli3(-exp(-μ / temp)) + PhaseGN.reli3(-exp(μ / temp))) / π
    if norm
        return (temp_independent + temp_dependent - vac_term) / temp^3
    else
        return temp_dependent + temp_independent - vac_term
    end
end

function der_omega(m, temp, mu, kappa, lambda)
    m = abs(m)
    temp_dependent = temp * m * (log(1 + exp((-m + mu) / temp)) + log(1 + exp((-m - mu) / temp))) / π
    temp_independent = (-m + m^2 - m^3 / (2 * lambda)) / π
    return temp_dependent + temp_independent - kappa
end

function find_mass(T, mu, kappa, lambda)
    f(x) = der_omega(x, T, mu, kappa, lambda)
    roots = find_zeros(f, (-2.0, 2.0))
    # omega.(roots, T, mu, kappa, lambda) 
    if omega(maximum(roots), T, mu, kappa, lambda) < 0
        return maximum(roots)
    else
        return 0.0
    end
end

function find_line(T, mu, lambda)
    f(x) = der_omega(x, T, mu, 0.0, lambda)
    roots = find_zeros(f, (-2.0, 2.0))
    if length(roots) == 1
        return 0.0
    else
        return 1.0
    end
end

function find_other_line(T, lambda)
    mu1 = find_critical_line(T, lambda)
    PhaseGN.find_zero_interval_start(x -> find_line(T, x, lambda), mu1, 1.5)
end

find_other_line(0.01, 5.0)
find_critical_line(0.01, 5.0)


function find_critical_line(T, lambda)
    mass(x) = find_mass(T, x, 0.0, lambda)
    PhaseGN.find_zero_interval_start(mass, 0.0, 1.4)
end

trange = sort(union(0.01:0.01:0.6, 0.61:0.005:0.7, 0.705:0.001:0.727))

mulist = find_critical_line.(trange, 5.0)

plot(mulist, trange)
plot!(crit_line.(trange), trange)
scatter!([crit_line(tcp(5.0))], [tcp(5.0)])
tcps = tcp(5.0)
trange2 = 0.01:0.01:tcp(5.0)
second_line = find_other_line.(trange2, 5.0)

plot!(second_line, trange2)

find_mass(0.1, 1.04, 0.0, 5.0)

p3 = @pgf Axis(
    {
        xlabel = L"\mu/M",
        ylabel = L"T/M",
        xmin = 0.00,
        ymin = 0.0,
        # ymax=0.79
    },
    Plot(
        {
        },
        Table(mulist, trange)
    ),
    Plot(
        {
            only_marks,
            mark_size = 1.5
        },
        Table([crit_line(tcp(5.0))], [tcp(5.0)])
    ),
    Plot(
        {
            style = "dashed"
        },
        Table(crit_line.(trange2), trange2)
    ),
    Plot(
        {
            style = "dashdotted"
        },
        Table(second_line, trange2)
    ),
)

pgfsave("phase_diagram_regulator.pdf", p3)

writedlm("fig_12_a_phase_diagram_cutoff.csv", hcat(trange2, crit_line.(trange2), second_line), ',')

writedlm("fig_12_a_critical_line_cutoff.csv", hcat(trange, mulist), ',')

function find_mass2(T, mu, kappa, lambda)
    f(x) = der_omega(x, T, mu, kappa, lambda)
    @show roots = find_zeros(f, (0.0, 2.0))
    omes = omega.(roots, T, mu, kappa, lambda)
    return roots[argmin(omes)]
end

function find_crit_line_kappa(T, kappa, lambda)
    function ff(μ)
        f(x) = der_omega(x, T, μ, kappa, lambda)
        roots = find_zeros(f, (0.0, 2.0))
        if length(roots) == 1
            return -1.0
        else
            return 1.0
        end
    end
    # @show ff(1.0)
    # murange = 0.9:0.001:1.3
    # data = zeros(length(murange))
    # for i in eachindex(murange)
    #     data[i]= ff(murange[i])
    # end
    # plot(murange, data, marker=:circle)
    @show fz = find_zeros(ff, 0.7, 1.5)
    if length(fz) == 2
        function ff2(μ)
            f(x) = der_omega(x, T, μ, kappa, lambda)
            roots = find_zeros(f, (0.0, 2.0))
            omes = omega.(roots, T, μ, kappa, lambda)
            return omes[1] - omes[end]
        end
        return find_zero(ff2, (fz[1], fz[2]))
    else
        return NaN
    end
end

find_crit_line_kappa(0.06, 0.046, 5.0)

trange3 = 0.01:0.001:0.14
mu_list2 = find_crit_line_kappa.(trange3, 0.01, 5.0)
plot!(mu_list2, trange3)

find_crit_line_kappa(0.066, 0.046, 5.0)

plot(murange, find_mass2.(0.11, murange, 0.01, 5.0), marker=:circle)

function find_pseudo_critical_line(mu, kappa, lambda)
    f(x) = find_mass2(x, mu, kappa, lambda)
    d2f(x) = sg_filter(f, x, dx=1e-4)
    plot(0.07:0.01:0.8, d2f, ylim=(-0.1, 1))
end

find_pseudo_critical_line(0.2, 0.0, 5.0)

mass_k(0.1, 1.0, 0.0, Parameters())
σ1(0.1, 1.0, Parameters())
find_mass2(0.1, 1.0, 0.01, 100.0)

find_zero(x -> der_omega(x, 0.1, 1.0, 0.01, 100.0), (0, 1.4))

Parameters()