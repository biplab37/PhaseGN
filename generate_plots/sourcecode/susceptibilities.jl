using Plots, LaTeXStrings, PhaseGN, PGFPlotsX

param2 = Parameters(Λ=5.0, κ=0.1)
param = Parameters(Λ=5.0, κ=kappa01(5.0))
kappa01(5.0)
M01 = (1 + sqrt(1 + 4 * π * 0.1)) / 2
function dm_dT(T, mu, param; dx=5e-2)
    return first_derivative(x -> mass_k(x, mu, 0.0, param), T, dx=dx)
end

function dm_dmu(T, mu, param; dx=1e-2)
    return first_derivative(x -> σ1(T, x, param), mu, dx=dx)
end

function d2m_dT2(T, mu, param; dx=5e-2)
    return sg_filter(x -> mass_k(x, mu, 0.0, param), T, dx=dx)
end

function d2m_dmu2(T, mu, param; dx=5e-2)
    return sg_filter(x -> mass_k(T, x, 0.0, param), mu, dx=dx)
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
mus = pseudo_critical_mu.(trange, param)

plot(trange, mus)

## determinant

function d2m_dt_dmu(temp, mu, param)
    f(t) = first_derivative(x -> mass_k(t, x, 0.0, param), mu)
    return first_derivative(f, temp)
end

function d2m_dmu_dt(temp, mu, param)
    f(m) = first_derivative(x -> mass_k(x, m, 0.0, param), temp)
    return first_derivative(f, mu)
end

function determinant(temp, mu, param)
    return d2m_dmu2(temp, mu, param) * d2m_dT2(temp, mu, param) - d2m_dmu_dt(temp, mu, param)^2
end

T1 = 0.1
murange = 0.0:0.01:1.2

function find_pseudo_critical(t, param)
    f(mu) = determinant(t, mu, param)
    return PhaseGN.bisection(f, 0.5, 1.13)
end

function find_pseudo_critical_mu(mu, param)
    f(t) = determinant(t, mu, param)
    return PhaseGN.bisection(f, 0.6, 0.75)
end
find_pseudo_critical(0.2, param)

Trange = 0.2:0.01:0.81

data_pseudo_critical = zeros(length(Trange))

Threads.@threads for i in eachindex(Trange)
    data_pseudo_critical[i] = find_pseudo_critical(Trange[i], param)
end

scatter(Trange, data_pseudo_critical, msw=0, ms=3)

scatter!([Trange[1]], [data_pseudo_critical[1]])

plot(murange, determinant.(0.75, murange, param))
vline!([M01])
hline!([0.0])
vline!(ans)


t1range = 0.01:0.01:0.1
t2range = 0.11:0.01:0.45
t3range = 0.46:0.01:0.7
t4range = 0.71:0.005:0.749
t5range = 0.749:0.005:0.82
t6range = 0.625:0.005:0.7

mu5range = 0.0:0.01:0.6
mu6range = 0.65:0.005:0.75

mu2 = find_pseudo_critical.(t2range, param)
scatter!(t2range, mu2)

mu3 = find_pseudo_critical.(t3range, param)
scatter!(t3range, mu3)

mu4 = find_pseudo_critical.(t4range, param)
scatter(t4range, mu4)

mu5 = find_pseudo_critical.(t5range, param)
scatter!(t5range, mu5)

t6 = all_zeros.(mu6range, param)
scatter!(t6, mu6range)

t5 = pseudo_critical_temp.(mu5range, 0.0, param)
scatter(t5, mu5range)

mu1 = find_pseudo_critical.(t2range, param)
scatter(t2range, mu2)

plot(determinant.(0.65, 0.6:0.001:0.9, param))

hline!([0.0])

temps = vcat(t2range, t3range, t4range, t5range)
mus = vcat(mu2, mu3, mu4, mu5)

ids = findall(x -> (abs(x - 0.65) < 1e-2), temps)

popat!(temps, 56)
popat!(mus, 56)

plot(temps, mus, marker=:circle)
scatter!(temps[55:56], mus[55:56])

p3 = @pgf Axis(
    {
}
)