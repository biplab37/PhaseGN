using PhaseGN, PGFPlotsX, LaTeXStrings, DelimitedFiles

data = readdlm("pressure_fluctuations.dat", skipstart=1)

trange = data[:, 1]
press = data[:, 2]

plot(trange, press, marker=:circle, xlabel=L"T/M", ylabel=L"P/T^3")

function phase_zero(ω, T, mu, param)
    impi = imagpart_phi_q(ω, T, mu, 0.0, param)
    repi = Π0_phi(T, mu, param) - PhaseGN.realpart(imagpart_phi_q, ω, T, mu, 0.0, param)
    return atan(impi, repi)
end
param = Parameters()

pressure_boost = zeros(length(trange))
Threads.@threads for i in eachindex(trange)
    pressure_boost[i] = pressure(phase_zero, trange[i], 0.0, param)
end

## pole approximation
function pole_approx(ω, T, mu, param)
    Mphi = M_phi(T, mu, param)
    m = σ1(T, mu, param)
    return (Mphi < ω < 2m) ? π : 0.0
end

pressure_pole = zeros(length(trange))

Threads.@threads for i in eachindex(trange)
    pressure_pole[i] = pressure(pole_approx, trange[i], 0.0, param)
end

function pole_approx_2(T, param)
    Mphi = M_phi(T, 0.0, param)
    integrand(q) = -q * log(1 - exp(-sqrt(q^2 + Mphi^2) / T)) / (2π * T^2)
    return PhaseGN.integrate(integrand, 0.0, param.Λ)
end

pressure_pole_2 = zeros(length(trange))

Threads.@threads for i in eachindex(trange)
    pressure_pole_2[i] = pole_approx_2(trange[i], param)
end

plot!(trange, pressure_boost, label="boost")
plot!(trange, pressure_pole, label="Pole approx")
plot!(trange, 2 * pressure_pole_2, label="Pole approx 2")

p_pole = @pgf Axis(
    {
        xlabel = L"T/M",
        ylabel = L"P/T^3",
        xmin = 0.0,
        xmax = 1.6,
        ymin = 0.0,
        # ymax = 0.35,
    },
    PlotInc(
        {
            color = "black",
            no_marks,
        },
        Table(x=trange, y=press)
    ),
    LegendEntry("Full"),
    PlotInc(
        {
            color = "black",
            no_marks,
            style = "dashed",
        },
        Table(x=trange, y=pressure_boost)
    ),
    LegendEntry("Boost"),
    PlotInc(
        {
            color = "black",
            no_marks,
            style = "dotted",
        },
        Table(x=trange, y=2 * pressure_pole_2)
    ),
    LegendEntry("Pole"),
)

pgfsave("pressure_comparison.pdf", p_pole)