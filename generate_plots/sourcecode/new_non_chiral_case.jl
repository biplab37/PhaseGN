using Plots, LaTeXStrings, PGFPlotsX, PhaseGN, Roots

# continuum limit Λ->∞

param = Parameters(κ=0.046, Λ=5.0)

function second_derivative(f, x; dx=1e-2)
    return (f(x + dx) + f(x - dx) - 2 * f(x)) / dx^2
end
function sg_filter(f, x; dx=3e-2)
    return (-1 * f(x - 2 * dx) + 16 * f(x - dx) - 30 * f(x) + 16 * f(x + dx) - 1 * f(x + 2 * dx)) / (12 * dx^2)
end
function first_derivative(f, x; dx=1e-2)
    return (-2 * f(x - 2 * dx) - f(x + dx) + f(x + dx) + 2 * f(x + 2 * dx)) / (10 * dx)
end

function dm_dt(T, mu, param)
    return first_derivative(x -> σ1(x, mu, param), T)
end

function d2m_dt2(T, mu, param, dx=5e-2)
    return sg_filter(x -> σ1(x, mu, param), T, dx=dx)
end

function plot_dm_dt(mu, param)
    trange = 0.05:0.01:0.8
    dms = dm_dt.(trange, mu, param)
    plot(trange, dms)
end

plot_dm_dt(1.12, param)

function plot_d2m_dt2(T, param)
    murange = 0.0:0.01:1.1
    dms = d2m_dt2.(T, murange, param)
    plot(murange, dms)
    hline!([0.0])
end

plot_d2m_dt2(0.77, param)

function find_pseudo_crit_mu(T, param, dx=5e-2)
    f(x) = d2m_dt2.(T, x, param, dx)
    return PhaseGN.bisection(f, 0.0, 1.1281)
end

trange = union(0.2:0.01:0.7, 0.71:0.001:0.8)
mus = find_pseudo_crit_mu.(trange, param)

trange2 = 0.05:0.01:0.19
mus2 = find_pseudo_crit_mu.(trange2, param, 1e-2)

plot(mus, trange)
plot!(mus2, trange2, marker=:circle)
scatter!([1.1281], [0.0])

crit_l = PhaseGN.critical_line.(0.03:0.005:0.73, param)

plot!(crit_l, 0.03:0.005:0.73)

tlist = union(0.0, trange2, trange)
mulist = union(1.1281, mus2, mus)

plot(mulist, tlist)

p1 = @pgf Axis(
    {
        xlabel = L"\mu/M",
        ylabel = L"T/M",
        xmin = 0.0,
        ymin = 0.0,
        legend_pos = "south west"
    },
    Plot(
        {
            style = "dashed"
        },
        Table(mulist, tlist)
    ),
    Plot(
        {
            style = "solid"
        },
        Table([1.0; crit_l], [0.0; collect(0.03:0.005:0.73)])
    ),
    Legend(
        [
        "Pseudocritical line",
        "Critical Line"
    ]
    )
)