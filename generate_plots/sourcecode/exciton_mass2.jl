function mphi(temp, mu, param)
    func(x) = PhaseGN.fullrealpart_phi(x, temp, mu, param)
    return find_zeros(func, 0.1, 3.0)
end

trange = 0.01:0.01:1.1

mass = σ1.(trange, 0.0, param)

Mphis2 = []
for i in 1:length(trange)
    mphis = mphi(trange[i], 0.0, param)
    push!(Mphis2, mphis)
end

plot(trange, [minimum.(Mphis2) maximum.(Mphis2)])


function delta_phi_at_threshold(q, T, param)
    m = mass_k(T, 0.0, q, param)
    thre = sqrt(q^2 + 4 * m^2) - min(1e-4, (sqrt(q^2 + 4 * m^2) - q) / 2)
    @show impi = (PhaseGN.imagpart_phi_q_refactored_m(thre, T, 0.0, q, m, param) + PhaseGN.imagpart_phi_q_refactored_m(thre, T, 0.0, q, m, param)) / 2
    @show repi = PhaseGN.realpart_3(PhaseGN.imagpart_phi_q_refactored_m, thre, T, 0.0, q, m, param)
    return atan(impi, repi)
end

delta_phi_at_threshold(0.0, 0.99, param)

function mmott_momenta2(T, param)
    func(x) = π - delta_phi_at_threshold(x, T, param)
    return PhaseGN.find_zero_interval_start(func, 0.0, 3.0)
end

mmott_momenta2(0.99, param)

function fullrealpart_phi2(ω, temp, μ, param)
    m = σ1(temp, μ, param)
    Ep(p) = sqrt(p^2 + m^2)
    integrand(p) = 1 / π + p * (1) * (1 - PhaseGN.numberF(temp, -μ, Ep(p)) - PhaseGN.numberF(temp, μ, Ep(p))) * (PhaseGN.PrincipalValue(ω - 2 * Ep(p)) - PhaseGN.PrincipalValue(ω + 2 * Ep(p))) / π
    return -1 / π + PhaseGN.integrate(integrand, Complex(0.0), Complex(param.Λ))
end

fullrealpart_phi(2 + 1im, 0.1, 0.0, param)

using UsefulFunctions
function reI2(omega, temp, mu, param)
    m = mass_k(temp, mu, 0.0, param)
    En(p) = sqrt(p^2 + m^2)
    integrand(p) = p * (1 - numberF(temp, mu, En(p)) - numberF(temp, -mu, En(p))) * (PrincipalValue(2 * En(p) - omega) + PrincipalValue(2 * En(p) + omega)) / (8π * En(p)^2)
    return PhaseGN.integrate(integrand, 0.0, param.Λ)
end

reI2(0.1, 0.01, 0.0, param)
function imI2(omega, temp, mu, param)
    m = mass_k(temp, mu, 0.0, param)
    if omega^2 <= 4 * m^2 || omega^2 >= 4 * (param.Λ^2 + m^2)
        return 0.0
    end
    return (1 - numberF(temp, mu, omega / 2) - numberF(temp, -mu, omega / 2)) / (8 * omega)
end

function I1(temp, mu, param)
    m = mass_k(temp, mu, 0.0, param)
    En(p) = sqrt(p^2 + m^2)
    integrand(p) = p * (1 - numberF(temp, mu, En(p)) - numberF(temp, -mu, En(p))) / (4π * En(p))
    return (param.Λ - 1.0) / (2π) - 2 * PhaseGN.integrate(integrand, 0.0, param.Λ)
end

I1(0.1, 0.0, param)
using NLsolve

function solve_two(temp, mu, param)
    m = mass_k(temp, mu, 0.0, param)
    i1 = I1(temp, mu, param)
    function func(x)
        rep = reI2(x[1], temp, mu, param)
        imp = imI2(x[1], temp, mu, param)
        term1 = x[1]^2 - 0.25 * x[2]^2 - 4 * m^2 - i1 * rep / (rep^2 + imp^2)
        term2 = x[1] * x[2]^2 - i1 * imp / (rep^2 + imp^2)
        return [term1, term2]
    end
    return nlsolve(func, [0.5, 0.0])
end

solve_two(0.9, 0.0, param)

function plot_two(temp, mu, param)
    m = mass_k(temp, mu, 0.0, param)
    i1 = I1(temp, mu, param)
    function func(x)
        rep = reI2(x[1], temp, mu, param)
        imp = imI2(x[1], temp, mu, param)
        term1 = x[1]^2 - 0.25 * x[2]^2 - 4 * m^2 - i1 * rep / (rep^2 + imp^2)
        term2 = x[1] * x[2]^2 - i1 * imp / (rep^2 + imp^2)
        return -term1 * term2
    end
    orange = 0.0:0.05:2.0
    grange = 0.0:0.01:2.0
    heatmap(orange, grange, (x, y) -> func([x, y]))
end

using Plots, PhaseGN
param = Parameters(Λ=5.0, κ=0.046)
plot_two(1.1, 0.0, param)
using Roots
function solve_two2(temp, mu, param)
    m = mass_k(temp, mu, 0.0, param)
    i1 = I1(temp, mu, param)

    rep(omega) = reI2(omega, temp, mu, param)
    imp(omega) = imI2(omega, temp, mu, param)
    gamma(omega) = i1 * imp(omega) / (omega * (rep(omega)^2 + imp(omega)^2))
    term1(omega) = omega^2 - 4 * m^2 - 0.25 * gamma(omega)^2 - i1 * rep(omega) / (rep(omega)^2 + imp(omega)^2)

    om1 = PhaseGN.bisection(term1, 1.0, 3.0)
    return om1, gamma(om1)
end
function solve_two2_phi(temp, mu, param)
    m = mass_k(temp, mu, 0.0, param)
    i1 = I1(temp, mu, param)

    rep(omega) = reI2(omega, temp, mu, param)
    imp(omega) = imI2(omega, temp, mu, param)
    gamma(omega) = i1 * imp(omega) / (omega * (rep(omega)^2 + imp(omega)^2))
    term1(omega) = omega^2 - 0.25 * gamma(omega)^2 - i1 * rep(omega) / (rep(omega)^2 + imp(omega)^2)

    om1 = PhaseGN.bisection(term1, 0.1, 1.8)
    return om1, gamma(om1)
end
solve_two2(1.0, 0.0, param)

trange5 = 0.01:0.01:1.3

new_mphis = zeros(length(trange5))
new_width_phi = zeros(length(trange5))
new_msigmas = zeros(length(trange5))
new_width_sigma = zeros(length(trange5))
Threads.@threads for i in eachindex(trange5)
    new_mphis[i], new_width_phi[i] = solve_two2_phi(trange5[i], 0.0, param)
    new_msigmas[i], new_width_sigma[i] = solve_two2(trange5[i], 0.0, param)
end

Mphis2 = []
for i in 1:length(trange5)
    mphis = mphi(trange5[i], 0.0, param)
    push!(Mphis2, mphis)
end

length(trange5)
mass = mass_k.(trange5, 0.0, 0.0, param)
σ1(0.01, 0.0, param)
mass_k(0.01, 0.0, 0.0, param)
plot(trange5, [2mass new_mphis new_msigmas minimum.(Mphis2) maximum.(Mphis2)], marker=:circle)

plot(trange, [2mass minimum.(Mphis2) maximum.(Mphis2)])
function reparti2(omega, gamma, temp, param)
    m = mass_k(temp, 0.0, 0.0, param)
    En(p) = sqrt(p^2 + m^2)
    integrand(p) = p * (1 - 2 * numberF(temp, 0.0, En(p))) * (1) * 2 * En(p) * (omega^2 - gamma^2 / 4 - En(p)^2) / ((omega^2 - gamma^2 / 4 - En(p)^2)^2 + gamma^2 * omega^2 / 4)
    return ((param.Λ - 1) + PhaseGN.integrate(integrand, 0.0, param.Λ)) / π
end

function imagparti2(omega, gamma, temp, param)
    m = mass_k(temp, 0.0, 0.0, param)
    En(p) = sqrt(p^2 + m^2)
    integrand(p) = p * (1 - 2 * numberF(temp, 0.0, En(p))) * (1) * 2 * En(p) * (omega * gamma) / ((omega^2 - gamma^2 / 4 - En(p)^2)^2 + gamma^2 * omega^2 / 4)
    return (PhaseGN.integrate(integrand, 0.0, param.Λ)) / π
end

function solve_two3(temp, param)

    function gamma(omega)
        rep(g) = reparti2(omega, g, temp, param)
        return find_zero(rep, -0.1)
    end
    f(omega) = imagparti2(omega, gamma(omega), temp, param)
    x = find_zero(f, 1.0)
    return x
end

solve_two2_phi(0.95, 0.0, param)
m = mass_k.(trange5, 0.0, 0.0, param)
using PGFPlotsX, LaTeXStrings
p = @pgf Axis(
    {
        xlabel = "T/M",
        ylabel = "Mass [M]",
        xmin = 0.0,
        xmax = 1.2,
        ymin = 0.0,
        minor_y_tick_num = 9,},
    PlotInc(
        {
            no_marks,
            black,
            dashed,
        },
        Table(trange5, 2 .* mass)
    ),
    LegendEntry("2m"),
    PlotInc(
        {
            no_marks,
            black,
        },
        Table(trange5, new_msigmas)
    ),
    LegendEntry(L"M_\sigma"),
    PlotInc(
        {
            no_marks,
            black,
            dashdotted,
        },
        Table(trange5, new_mphis)
    ),
    LegendEntry(L"M_\varphi"),
)

writedlm("fig_4_a_exciton_masses.csv", hcat(trange5, 2 .*mass, new_msigmas, new_mphis), ',')

pgfsave("exciton_masses.pdf", p)

plot(trange5, [new_width_phi new_width_sigma])

p2 = @pgf Axis(
    {
        xlabel = L"T/M",
        ylabel = "Width [M]",
        xmin = 0.0,
        xmax = 1.2,
        legend_pos = "north west"
    },
    PlotInc(
        {
            no_marks,
            color = "black",
            style = "dashdotted",
        },
        Table(trange5, new_width_phi)
    ),
    LegendEntry(L"\Gamma_\varphi"),
    PlotInc(
        {
            no_marks,
            color = "black",
            style = "solid",
        },
        Table(trange5, new_width_sigma),
    ),
    LegendEntry(L"\Gamma_\sigma"))

pgfsave("width.pdf", p2)

writedlm("fig_4_b_widths.csv", hcat(trange5, new_width_phi, new_width_sigma), ',')

delta_phi_at_threshold(0.0, 0.999, param)