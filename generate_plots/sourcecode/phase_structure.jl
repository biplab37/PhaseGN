using PhaseGN, Plots, LaTeXStrings, DelimitedFiles

kappa01(5.0)
param = Parameters(Λ=5.0, κ=0.01)
function plot_omega_kappa(T, mu, param)
    σrange = -2.0:0.01:2.0
    omega_kappa = [PhaseGN.Omega_kappa(σ, T, mu, param) for σ in σrange]
    plot(σrange, omega_kappa, label="T=$T, μ=$mu", xlabel=L"\sigma", ylabel=L"\Omega(\sigma, T, \mu)", legend=:topleft)
end

plot_omega_kappa(0.01, 1.0, param)

Trange = 0.01:0.01:1.3
murange = 0.0:0.01:1.3

phase_diag = zeros(length(Trange), length(murange))

Threads.@threads for i in eachindex(Trange)
    for j in eachindex(murange)
        phase_diag[i, j] = PhaseGN.mass_k(Trange[i], murange[j], 0.0, param)
    end
end
plotly()
surface(Trange, murange, phase_diag, xlabel=L"T", ylabel=L"\mu", title="Phase diagram", legend=:topleft)

σrange = -2.0:0.01:2.0

omega_phase_t = zeros(length(σrange), length(murange))

param2 = Parameters(Λ=5.0, κ=0.01)
Threads.@threads for i in eachindex(murange)
    for j in eachindex(σrange)
        omega_phase_t[j, i] = PhaseGN.Omega_kappa(σrange[j], 0.005, murange[i], param)
    end
end

surface(σrange, murange, omega_phase_t', xlabel=L"\sigma", ylabel=L"\Omega(\sigma, T, \mu)", title="Phase diagram", legend=:topleft)

plot_phase_diag(param)