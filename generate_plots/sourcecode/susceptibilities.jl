using Plots, LaTeXStrings, PhaseGN, PGFPlotsX

param2 = Parameters(Λ=5.0, κ=0.1)
param = Parameters(Λ=5.0, κ=0.046)
M01 = (1 + sqrt(1 +4*π*0.1))/2
function dm_dT(T, mu, param; dx=5e-2)
    return first_derivative(x->mass_k(x, mu, 0.0, param), T, dx=dx)
end

function dm_dmu(T, mu, param; dx=1e-2)
    return first_derivative(x->σ1(T, x, param), mu, dx=dx)
end

function d2m_dT2(T, mu, param; dx=5e-2)
    return sg_filter(x->mass_k(x, mu, 0.0, param), T, dx=dx)
end

function d2m_dmu2(T, mu, param; dx=5e-2)
    return sg_filter(x->mass_k(T, x, 0.0, param), mu, dx=dx)
end

function plot_dm_dmu(temp, param)
    murange = 0.0:0.01:1.5
    dm = [dm_dmu(temp, mu, param) for mu in murange]
    plot(murange, dm)
    vline!([M01])
end

function plot_dm_dT(temp, param)
    Trange = 0.1:0.01:1.5
    dm = [dm_dT(temp, mu, param) for temp in Trange]
    plot(Trange, dm)
end

function plot_mass_mu(temp, param)
    murange = 0.0:0.01:1.5
    masses = mass_k.(temp, murange, 0.0, param)
    plot(murange, masses)
end

plot_dm_dmu(0.9, param)
plot_mass_mu(1.0, param)

function pseudo_critical_mu(temp, param)
    der2_mass(mu, param) = sg_filter(x -> mass_k(temp, x, 0.0, param), mu, dx=5e-2)
    try
        return find_zero(x -> der2_mass(x, param), (0.01, 1.5))
    catch
        return 0.0
    end
end

trange = 0.1:0.01:1.0
mus = pseudo_critical_mu.(trange, param2)

plot(trange, mus)

## determinant

function d2m_dt_dmu(temp, mu, param)
    f(t) = first_derivative(x->mass_k(t, x, 0.0, param), mu)
    return first_derivative(f, temp)
end

function d2m_dmu_dt(temp, mu, param)
    f(m) = first_derivative(x->mass_k(x, m, 0.0, param), temp)
    return first_derivative(f, mu)
end

function determinant(temp, mu, param)
    return d2m_dmu2(temp, mu, param)*d2m_dT2(temp, mu, param) - d2m_dmu_dt(temp, mu, param)^2
end

T1 = 0.1
murange = 0.:0.01:1.2

function find_pseudo_critical(t, param)
    f(mu) = determinant(t, mu, param)
    return PhaseGN.bisection(f, 0.0, 1.3)
end

find_pseudo_critical(0.2, param)

Trange = 0.2:0.01:0.9

data_pseudo_critical = zeros(length(Trange))

Threads.@threads for i in eachindex(Trange)
    data_pseudo_critical[i] = find_pseudo_critical(Trange[i], param2)
end

scatter!(Trange, data_pseudo_critical, msw=0, ms=3)

scatter!([Trange[1]], [data_pseudo_critical[1]])

plot(murange, determinant.(0.2, murange, param))
vline!([M01])
hline!([0.0])
vline!(ans)

