## This file contains phase shift definitions with the refactored version of code used in the paper.

function phase_shift(imagpart::Function, Π0::Function, mass::Function, ω, q, T, mu, param::Parameters)
    m = mass(T, mu, q, param)
    impi = imagpart(ω, T, mu, q, m, param)
    repi = Π0(T, mu, param) - realpart(imagpart, ω, T, mu, q, m, param)
    return atan(impi, repi)
end

phase_shift_phi(ω, q, T, mu, param) = phase_shift(imagpart_phi_q_refactored_m, Π0_phi, mass_k, ω, q, T, mu, param)
phase_shift_sigma(ω, q, T, mu, param) = phase_shift(imagpart_sigma_q_refactored_m, Π0_sigma, mass_k, ω, q, T, mu, param)