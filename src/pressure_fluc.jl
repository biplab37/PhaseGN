## Fluctuation pressure due to non-zero external momentum

function pressure_fluc(phase::Function, temp, μ, param::Parameters)
    phi(ω, q) = phase(temp, μ, ω, q, param)
    integrand(ω, q) = q * (1.0 + 2.0 / (exp(ω / temp) - 1.0)) * phi(ω, q) / (4π^2)
    return integrate(x -> integrand(x[1], x[2]), [0.0, 0.0], [2 * param.Λ, param.Λ])
end

function pressure_fluc_phi(temp, μ, param::Parameters)
    return pressure_fluc(phasetot_phi, temp, μ, param)
end

function pressure_fluc_sigma(temp, μ, param::Parameters)
    return pressure_fluc(phasetot_sigma, temp, μ, param)
end

export pressure_fluc_sigma, pressure_fluc_phi, pressure_fluc
