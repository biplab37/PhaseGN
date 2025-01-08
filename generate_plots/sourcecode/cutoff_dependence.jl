using Plots, PhaseGN, PGFPlotsX, LaTeXStrings, DelimitedFiles

param = Parameters(Λ=2.0)

function delta_phi(ω, q, T, param)
    impi = imagpart_phi_q(ω, T, 0.0, q, param)
    repi = Π0_phi(T, 0.0, param) - PhaseGN.realpart(imagpart_phi_q, ω, T, 0.0, q, param)
    return atan(impi, repi)
end

function phase_shift_data(ωrange, q, T, param)
    phase_shift = delta_phi.(ωrange, q, T, param)
    return phase_shift
end

# plot the phase shift
q = 1.0
T = 0.85
ωrange = sort(union(10 .^ (-2:0.05:1.2), [sqrt(4 * σ1(T, 0.0, param)^2 + q^2)]))

phase_shifts = phase_shift_data(ωrange, q, T, param)

plot(ωrange, phase_shifts, xaxis=:log, xlabel=L"\omega", ylabel=L"\phi_\varphi", label="q=$q")