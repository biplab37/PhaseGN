## Entropy density calculations from the Beth-Uhlenbeck and the generalized Beth-Uhlenbeck formula.
function generalized_phase_shift(phase)
    return phase - sin(2 * phase) / 2
end

function entropy_density_fluc_phi_BU(T, mu, param; LD_factor=1.0)
    delta(om, q) = phase_shift_phi(om, q, T, mu, param)
    sigma(om, T) = om / (4 * T^2 * sinh(om / 2T)^2)
    integrand(om, q) = q * sigma(om, T) * delta(om, q) / (4 * π^2)

    cutoff = 2 * param.Λ * sqrt(2 + LD_factor^2 / 4)

    return integrate(x -> integrand(x...), [0.0, 0.0], [cutoff, LD_factor * param.Λ])
end

function entropy_density_fluc_phi_genBU(T, mu, param; LD_factor=1.0)
    delta(om, q) = generalized_phase_shift(phase_shift_phi(om, q, T, mu, param))
    sigma(om, T) = om / (4 * T^2 * sinh(om / 2T)^2)
    integrand(om, q) = q * sigma(om, T) * delta(om, q) / (4 * π^2)

    cutoff = 2 * param.Λ * sqrt(2 + LD_factor^2 / 4)

    return integrate(x -> integrand(x...), [0.0, 0.0], [cutoff, LD_factor * param.Λ])
end

#TODO: Update the formulas of the number densities in the functions below!!
function number_density_fluc_phi_BU(T, mu, param; LD_factor=1.0)
    delta(om, q) = phase_shift_phi(om, q, T, mu, param)
    sigma(om, T) = om / (4 * T^2 * sinh(om / 2T)^2)
    integrand(om, q) = q * sigma(om, T) * delta(om, q) / (4 * π^2)

    cutoff = 2 * param.Λ * sqrt(2 + LD_factor^2 / 4)

    return integrate(x -> integrand(x...), [0.0, 0.0], [cutoff, LD_factor * param.Λ])
end

function number_density_fluc_phi_genBU(T, mu, param; LD_factor=1.0)
    delta(om, q) = generalized_phase_shift(phase_shift_phi(om, q, T, mu, param))
    sigma(om, T) = om / (4 * T^2 * sinh(om / 2T)^2)
    integrand(om, q) = q * sigma(om, T) * delta(om, q) / (4 * π^2)

    cutoff = 2 * param.Λ * sqrt(2 + LD_factor^2 / 4)

    return integrate(x -> integrand(x...), [0.0, 0.0], [cutoff, LD_factor * param.Λ])
end

export entropy_density_fluc_phi_BU, entropy_density_fluc_phi_genBU
