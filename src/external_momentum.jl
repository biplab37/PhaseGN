function imagpart_sigma_q(ω, temp, μ, q, param)
    if q == 0.0
        return imagpart_sigma(ω, temp, μ, param)
    end
    m = σ1(temp, μ, param)
    s = ω^2 - q^2
    if s < 0.0
        y = q * sqrt(1 - 4m^2 / s)
        pauli_blocking3(E1) = numberF(temp, μ, 0.5 * (E1 - ω)) - numberF(temp, μ, 0.5 * (E1 + ω)) + numberF(temp, μ, 0.5 * (-E1 - ω)) - numberF(temp, μ, 0.5 * (-E1 + ω))
        integrand3(E1) = -(s - 4 * m^2) * pauli_blocking3(E1) / (sqrt(4 * ((E1^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E1^2 * ω^2) * 32π)
        return integrate(integrand3, y, param.Λ)
    end
    if s > 4m^2
        y = q * sqrt(1 - 4m^2 / s)
        if ω > 0.0
            pauli_blocking1(E2) = numberF(temp, μ, 0.5 * (-E2 - ω)) - numberF(temp, μ, 0.5 * (-E2 + ω))
            integrand1(E2) = (s - 4 * m^2) * pauli_blocking1(E2) / (sqrt(4 * ((E2^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E2^2 * ω^2) * 32π)
            return integrate(integrand1, -y, y)
        else
            pauli_blocking2(E2) = numberF(temp, μ, 0.5 * (E2 - ω)) - numberF(temp, μ, 0.5 * (E2 + ω))
            integrand2(E2) = (s - 4 * m^2) * pauli_blocking2(E2) / (sqrt(4 * ((E2^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E2^2 * ω^2) * 32π)
            return integrate(integrand2, -y, y)
        end
    end
    return 0.0
end

function imagpart_phi_q(ω, temp, μ, q, param)
    if q == 0.0
        return imagpart_phi(ω, temp, μ, param)
    end
    m = σ1(temp, μ, param)
    s = ω^2 - q^2
    if s < 0.0
        y = q * sqrt(1 - 4m^2 / s)
        pauli_blocking3(E1) = numberF(temp, μ, 0.5 * (E1 - ω)) - numberF(temp, μ, 0.5 * (E1 + ω)) + numberF(temp, μ, 0.5 * (-E1 - ω)) - numberF(temp, μ, 0.5 * (-E1 + ω))
        integrand3(E1) = -(s) * pauli_blocking3(E1) / (sqrt(4 * ((E1^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E1^2 * ω^2) * 32π)
        return integrate(integrand3, y, param.Λ)
    end
    if s > 4m^2
        y = q * sqrt(1 - 4m^2 / s)
        if ω > 0.0
            pauli_blocking1(E2) = numberF(temp, μ, 0.5 * (-E2 - ω)) - numberF(temp, μ, 0.5 * (-E2 + ω))
            integrand1(E2) = (s) * pauli_blocking1(E2) / (sqrt(4 * ((E2^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E2^2 * ω^2) * 32π)
            return integrate(integrand1, -y, y)
        else
            pauli_blocking2(E2) = numberF(temp, μ, 0.5 * (E2 - ω)) - numberF(temp, μ, 0.5 * (E2 + ω))
            integrand2(E2) = (s) * pauli_blocking2(E2) / (sqrt(4 * ((E2^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E2^2 * ω^2) * 32π)
            return integrate(integrand2, -y, y)
        end
    end
    return 0.0
end

function Π0_phi(temp, μ, q, param)
    m = σ1(temp, μ, param)
    integrand(Ep) = 1 / π - Ep * (1) * (1 - numberF(temp, -μ, Ep) - numberF(temp, μ, Ep)) / (π * Ep)
    return -1 / π + integrate(integrand, m, sqrt(param.Λ^2 + m^2))
end

function Π0_sigma(temp, μ, q, param)
    m = σ1(temp, μ, param)
    integrand(Ep) = 1 / π - Ep * (1 - m^2 / Ep^2) * (1 - numberF(temp, -μ, Ep) - numberF(temp, μ, Ep)) / (π * Ep)
    return -1 / π + integrate(integrand, m, sqrt(param.Λ^2 + m^2))
end

export imagpart_sigma_q, imagpart_phi_q, Π0_phi, Π0_sigma