using PhaseGN, PGFPlotsX, LaTeXStrings, DelimitedFiles

param = Parameters(Λ=2.0)

## This file contains code to explore the phase shift with momentum dependent mass gap

using Plots

function delta_phi(ω, q, T, param)
    impi = imagpart_phi_q(ω, T, 0.0, q, param)
    repi = Π0_phi(T, 0.0, param) - PhaseGN.realpart(imagpart_phi_q, ω, T, 0.0, q, param)
    return atan(impi, repi)
end

function delta_phi_refactored(ω, q, T, param)
    impi = PhaseGN.imagpart_phi_q_refactored_k(ω, T, 0.0, q, param)
    repi = Π0_phi(T, 0.0, param) - PhaseGN.realpart(PhaseGN.imagpart_phi_q_refactored_k, ω, T, 0.0, q, param)
    return atan(impi, repi)
end


function phase_shift_data(ωrange, q, T, param)
    phase_shift = delta_phi.(ωrange, q, T, param)
    return phase_shift
end

function phase_shift_data_refactored(ωrange, q, T, param)
    phase_shift = delta_phi_refactored.(ωrange, q, T, param)
    return phase_shift
end

q = 5.0
T = 0.85

param = Parameters()

m, mk = σ1(T, 0.0, param), mass_k(T, 0.0, q, param)

sqrt(4 * m^2 + q^2), sqrt(4 * mk^2 + q^2)

ωrange = (q-0.1):0.01:(q+0.2)

phase_shifts = phase_shift_data(ωrange, q, T, param)
phase_shifts_refactored = zeros(length(ωrange))
@time Threads.@threads for i in eachindex(ωrange)
    phase_shifts_refactored[i] = delta_phi_refactored(ωrange[i], q, T, param)
end

plot(ωrange, phase_shifts, xaxis=:log, xlabel=L"\omega", ylabel=L"\phi_\varphi", label="original")
plot!(ωrange, phase_shifts_refactored, xaxis=:log, xlabel=L"\omega", ylabel=L"\phi_\varphi", label="refactored")

q_list = [0.0, 1.0, 2.0]

threshold_values = [sqrt(4 * mass_k(T, 0.0, q, param)^2 + q^2) for q in q_list]

ωrange = sort(union(10 .^(-2:0.1:1.), threshold_values, [q-0.1:0.01:q+0.2 for q in q_list[2:end]]...))

phases_q = zeros(length(q_list), length(ωrange))
using ProgressMeter

@showprogress Threads.@threads for i in eachindex(ωrange)
    phases_q[1, i] = delta_phi(ωrange[i], q_list[1], T, param)
    phases_q[2, i] = delta_phi_refactored(ωrange[i], q_list[2], T, param)
    phases_q[3, i] = delta_phi_refactored(ωrange[i], q_list[3], T, param)
end

threshold_values2 = [sqrt(4 * σ1(T, 0.0, param)^2 + q^2) for q in q_list]
ωrange2 = sort(union(0.0:0.05:5.0, [threshold_values2[i]-0.1:0.005:threshold_values2[i]+0.1 for i in eachindex(q_list)]...))
phases_q0 = zeros(length(q_list), length(ωrange2))

@showprogress Threads.@threads for i in eachindex(ωrange2)
    phases_q0[1, i] = delta_phi(ωrange2[i], q_list[1], T, param)
    phases_q0[2, i] = delta_phi(ωrange2[i], q_list[2], T, param)
    phases_q0[3, i] = delta_phi(ωrange2[i], q_list[3], T, param)
end


plot(ωrange.^2 , phases_q[1, :], xlabel=L"s=\omega^2 - q^2", ylabel=L"\phi_\varphi", label="q=0.0")
plot!(ωrange.^2 .- q_list[2]^2 , phases_q[2, :], label="q=2.0")
plot!(ωrange.^2 .- q_list[3]^2 , phases_q[3, :], label="q=5.0", xlims=(0.0, 25.0))
savefig("phase_shift_with_mom_dep_mass_gap.pdf")


plot(ωrange2.^2 , phases_q0[1, :], xlabel=L"s=\omega^2 - q^2", ylabel=L"\phi_\varphi", label="q=0.0")
plot!(ωrange2.^2 .- q_list[2]^2 , phases_q0[2, :], label="q=2M")
plot!(ωrange2.^2 .- q_list[3]^2 , phases_q0[3, :], label="q=5M")
savefig("phase_shift_without_momentum_dependent_mass_gap_.pdf")