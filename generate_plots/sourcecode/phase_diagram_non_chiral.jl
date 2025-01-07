using Plots, PhaseGN, UsefulFunctions, LaTeXStrings
set_theme(:default)
scalefontsizes(1.5)

param = Parameters(Λ=5.0, κ=0.046)
param2= Parameters(Λ=5.0, κ=0.1)
param3= Parameters(Λ=5.0, κ=1.)

mu_list = 1.0:0.01:1.2
mu_list2 = 0.0:0.01:1.25
trange = 0.01:0.005:1.2
length(mu_list2)
using Roots
function second_derivative(f, x; dx=2e-2)
    return (f(x+dx) + f(x-dx) - 2*f(x))/dx^2
end
function sg_filter(f, x; dx=4e-2)
    return (-1*f(x-2*dx) + 16*f(x-dx) - 30*f(x) + 16*f(x+dx) - 1*f(x+2*dx))/(12*dx^2)
end
function first_derivative(f, x; dx=1e-2)
    return (-2*f(x-2*dx) - f(x+dx) + f(x+dx) + 2*f(x+2*dx))/(10*dx)
end

function pseudo_critical_temp(mu, q, param)
    der2_mass(t, param) = sg_filter(x->mass_k(x, mu, q, param), t)
    try
        return find_zero(x->der2_mass(x, param), (0.1,1.5))
    catch
        return  0.0
    end
end

masses = mass_k.(trange, 1.2, 0.0, param)
dmdT = [first_derivative(x->mass_k(x, 1.2, 0.0, param), t) for t in trange]
plot(trange, masses)
plot!(trange, dmdT)

m0 = 0.1*4/π

findmin(dmdT)

function plot_m_dm(mu)
    masses = mass_k.(trange, mu, 0.0, param2)
    dmdT = [first_derivative(x -> mass_k(x, mu, 0.0, param2), t) for t in trange]
    plot(trange, masses)
    plot!(trange, dmdT)
end

plot_m_dm(1.)

function find_pseudo_critical(mu, param; trange = trange)
    dmdT = [first_derivative(x -> mass_k(x, mu, 0.0, param), t) for t in trange]
    _, index = findmin(dmdT)
    return trange[index]
end

temps = find_pseudo_critical.(mu_list, param)

plot(temps, mu_list)

plot_m_dm(1.13)

temperatures = pseudo_critical_temp.(mu_list, 0.0, param2)
temperatures1 = pseudo_critical_temp.(mu_list, 0.0, param2)
temperatures2 = pseudo_critical_temp.(mu_list2, 0., param)
temperatures3 = pseudo_critical_temp.(mu_list2, .0, param3)

mu_range = PhaseGN.critical_line.(trange, param)

scalefontsizes(1.3)
plot(temperatures1, mu_list2, marker=:circle, xlabel=L"T/M", ylabel=L"\mu/M", label=L"\kappa=0.1 M^2", fontfamily="Computer Modern")
# plot!(temperatures2[1:end-13], mu_list2[1:end-13])
plot!(trange,  mu_range, lab=L"\kappa=0")
hline!([mu_c(param2)], opacity=0.5, l=:dash, color=:black, lab="")
plot!(xlims=(0.1, 0.95))
plot!(temperatures, mu_list)

savefig("pseudo_critical_phase_transition.pdf")
plot!(temperatures3, mu_list2)



function plot_omega(T, mu, param)
    mrange =  -1.5:0.01:1.5
    omegas = PhaseGN.Omega_kappa.(mrange, T, mu, param)
    plot(mrange, omegas) 
end

function plot_omega_t(trange, mu, param)
    mrange = 0.0:0.01:1.5
    omegas = zeros(length(trange), length(mrange))
    Threads.@threads for i in eachindex(trange)
        omegas[i, :] = PhaseGN.Omega_kappa.(mrange, trange[i], mu, param)
    end
    plot(mrange, omegas', labels=trange')
end

plot_omega_t(0.01:0.05:0.21, 1.13, param)

masses = mass_k.(0.01:0.01:1.0, 1.13, 0.0, param)
plot(masses)

plot_omega(0.9, .1, param2)

function mu_c(param)
    return param.M*(1 + sqrt(1 + 4*param.κ*π/param.M^2))/2
end

mu_c(param)


pseudo_critical_temp(0.0, 0.0, Parameters(Λ=5.0, κ=0.0))
#/mu_c(Parameters(Λ=5.0, κ=0.04))
1/(2*log(2))

mass_k(0.72135, 0.0, 0.0, Parameters(Λ=5.0, κ=0.0))

p = plot_m_dm(1.05)

mu_range1 = 0.9:0.01:1.25

data_dmdt = zeros(length(mu_range1), length(trange))
data_masses = zeros(length(mu_range1), length(trange))

Threads.@threads for i in eachindex(mu_range1)
    # data_dmdt[i, :] = [first_derivative(x -> mass_k(x, mu_range1[i], 0.0, param2), t) for t in trange]
    data_masses[i, :] = [mass_k(t, mu_range1[i], 0.0, param2) for t in trange]
end

anim = @animate for i in eachindex(mu_range1)
    plot(trange, data_dmdt[i,:], lab="μ = $(mu_range1[i])M", legend=:topright, xlabel="T/M", ylabel=L"\frac{\partial m}{\partial T}")
    val, index = findmin(data_dmdt[i,:])
    scatter!([trange[index]], [val], lab="minima")
end
gif(anim, "dmdT.gif",fps=3)

anim2 = @animate for i in eachindex(mu_range1)
    plot(trange, data_masses[i,:], lab="μ = $(mu_range1[i])M", legend=:topright, xlabel="T/M", ylabel=L"m/M", ylims=(0.2,1.3))
    # val, index = findmin(data_dmdt[i,:])
    # scatter!([trange[index]], [val], lab="minima")
end
gif(anim2, "masses.gif",fps=3)
