## Fluctuation pressure due to non-zero external momentum

function pressure_fluc(phase::Function, temp, μ, param::Parameters)
    phi(ω, q) = phase(temp, μ, ω, q, param)

    #TODO: Check the formula. Seems to be wrong
    integrand(ω, q) = q * (1.0 / (exp(ω / temp) - 1.0)) * phi(ω, q) / (4π^2)
    return integrate(x -> integrand(x[1], x[2]), [0.001, 0.001], [2.0, 10.0], maxevals=1000)
end

function pressure_fluc_phi(temp, μ, param::Parameters)
    return pressure_fluc(phasetot_phi_q, temp, μ, param)
end

function pressure_fluc_sigma(temp, μ, param::Parameters)
    return pressure_fluc(phasetot_sigma_q, temp, μ, param)
end

export pressure_fluc_sigma, pressure_fluc_phi, pressure_fluc, pressure_fluc_phi_fast
