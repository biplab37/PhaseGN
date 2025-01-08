using Plots, PhaseGN, PGFPlotsX, LaTeXStrings

mu1 = 1.13 # this gave 3 zeros for susceptibility

trange = 0.01:0.005:0.7

param2 = Parameters(Λ=5.0, κ=0.1)

PhaseGN.pressure_MF_kappa(0.1, 0.0, param2)

function number_density(temp, mu, param)
    P(t) = PhaseGN.pressure_MF_kappa(t, mu, param)
    return first_derivative(P, temp, dx=5e-3)
end

pres = [PhaseGN.pressure_MF_kappa(t, mu1, param2, norm=true) for t in trange]
num = [number_density(t, mu1, param2) / t^2 for t in trange]

plot(trange, [pres num])