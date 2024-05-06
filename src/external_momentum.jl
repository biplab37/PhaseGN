"""
    imagpart_sigma_q(ω, temp, μ, q, param)

Returns the imaginary part of the polarisation function of the scalar channel for a given frequency ω, temperature temp, chemical potential μ, and momentum q.
"""
function imagpart_sigma_q(ω, temp, μ, q, param)
    if q == 0.0
        return imagpart_sigma(ω, temp, μ, param)
    else
        m = σ1(temp, μ, param)
        s = ω^2 - q^2

        if abs(ω) > 2 * sqrt(param.Λ^2 + m^2)
            return 0.0

        elseif s < 0.0
            y = q * sqrt(1 - 4m^2 / s)
            pauli_blocking3(E1) = numberF(temp, μ, 0.5 * (E1 - ω)) - numberF(temp, μ, 0.5 * (E1 + ω)) + numberF(temp, μ, 0.5 * (-E1 - ω)) - numberF(temp, μ, 0.5 * (-E1 + ω))
            integrand3(E1) = -(s - 4 * m^2) * pauli_blocking3(E1) / (sqrt(4 * ((E1^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E1^2 * ω^2) * 4π)
            return integrate(integrand3, y, sqrt(param.Λ^2 + m^2))

        elseif s > 4m^2
            y = q * sqrt(1 - 4m^2 / s)

            if ω > 0.0
                pauli_blocking1(E2) = numberF(temp, μ, 0.5 * (-E2 - ω)) - numberF(temp, μ, 0.5 * (-E2 + ω))
                integrand1(E2) = (s - 4 * m^2) * pauli_blocking1(E2) / (sqrt(4 * ((E2^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E2^2 * ω^2) * 4π)
                return integrate(integrand1, -y, y)
            else
                pauli_blocking2(E2) = numberF(temp, μ, 0.5 * (E2 - ω)) - numberF(temp, μ, 0.5 * (E2 + ω))
                integrand2(E2) = (s - 4 * m^2) * pauli_blocking2(E2) / (sqrt(4 * ((E2^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E2^2 * ω^2) * 4π)
                return integrate(integrand2, -y, y)
            end

        else
            return 0.0
        end
    end
end

"""
    imagpart_phi_q(ω, temp, μ, q, param)

Returns the imaginary part of the polarisation function for the pseudo-scalar channel for a given frequency ω, temperature temp, chemical potential μ, and momentum q.
"""
function imagpart_phi_q(ω, temp, μ, q, param)
    if q == 0.0
        return imagpart_phi(ω, temp, μ, param)
    else
        m = σ1(temp, μ, param)
        s = ω^2 - q^2
        if abs(ω) >= 2 * sqrt(param.Λ^2 + m^2)
            return 0.0

        elseif s < 0.0
            y = q * sqrt(1 - 4m^2 / s)
            pauli_blocking3(E1) = numberF(temp, μ, 0.5 * (E1 - ω)) - numberF(temp, μ, 0.5 * (E1 + ω)) + numberF(temp, μ, 0.5 * (-E1 - ω)) - numberF(temp, μ, 0.5 * (-E1 + ω))
            integrand3(E1) = -(s) * pauli_blocking3(E1) / (sqrt(4 * ((E1^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E1^2 * ω^2) * 4π)
            return integrate(integrand3, y, sqrt(param.Λ^2 + m^2))

        elseif s > 4m^2
            y = q * sqrt(1 - 4m^2 / s)

            if ω > 0.0
                pauli_blocking1(E2) = numberF(temp, μ, 0.5 * (-E2 - ω)) - numberF(temp, μ, 0.5 * (-E2 + ω))
                integrand1(E2) = (s) * pauli_blocking1(E2) / (sqrt(4 * ((E2^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E2^2 * ω^2) * 4π)
                return integrate(integrand1, -y, y)
            else
                pauli_blocking2(E2) = numberF(temp, μ, 0.5 * (E2 - ω)) - numberF(temp, μ, 0.5 * (E2 + ω))
                integrand2(E2) = (s) * pauli_blocking2(E2) / (sqrt(4 * ((E2^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E2^2 * ω^2) * 4π)
                return integrate(integrand2, -y, y)
            end

        else
            return 0.0
        end
    end
end

Π0_phi(temp, μ, q, param) = Π0_phi(temp, μ, param)

Π0_sigma(temp, μ, q, param) = Π0_sigma(temp, μ, param)

export imagpart_sigma_q, imagpart_phi_q, Π0_phi, Π0_sigma
