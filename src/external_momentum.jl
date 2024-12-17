"""
    imagpart_sigma_q(ω, temp, μ, q, param)

Returns the imaginary part of the polarisation function of the scalar channel for a given frequency ω, temperature temp, chemical potential μ, and momentum q.
"""
function imagpart_sigma_q(ω, temp, μ, q, param)
    if q == 0.0
        return imagpart_sigma(ω, temp, μ, param)
    end

    m = σ1(temp, μ, param)
    s = ω^2 - q^2

    if abs(s) > 2 * sqrt(param.Λ^2 + m^2)
        return 0.0
    end

    if (0.0 <= s <= 4m^2)
        return 0.0
    end

    if s < 0.0
        y = q * sqrt(1 - 4m^2 / s)
        pauli_blocking3(E1) = numberF(temp, μ, 0.5 * (E1 - ω)) - numberF(temp, μ, 0.5 * (E1 + ω)) + numberF(temp, μ, 0.5 * (-E1 - ω)) - numberF(temp, μ, 0.5 * (-E1 + ω))
        integrand3(E1) = -(s - 4 * m^2) * pauli_blocking3(E1) / (sqrt(4 * ((E1^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E1^2 * ω^2) * 4π)
        return integrate(integrand3, y, sqrt(param.Λ^2 + m^2))

    else # s > 4m^2
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
    end
end

function imagpart_sigma_q_tinf(ω, temp, μ, q, param)
    if q == 0.0
        return imagpart_sigma(ω, temp, μ, param)
    else
        m = 0.0
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
function imagpart_phi_q1(ω, temp, μ, q, param)
    if q == 0.0
        return imagpart_phi(ω, temp, μ, param)
    else
        m = σ1(temp, μ, param)
        s = ω^2 - q^2
        if abs(ω) >= 2 * sqrt(param.Λ^2 + m^2)
            return 0.0

        elseif s < 0.0
            y = q * sqrt(1 - 4m^2 / s)
            if y >= sqrt(param.Λ^2 + m^2)
                return 0.0
            end
            pauli_blocking3(E1) = numberF(temp, μ, 0.5 * (E1 - ω)) - numberF(temp, μ, 0.5 * (E1 + ω)) + numberF(temp, μ, 0.5 * (-E1 - ω)) - numberF(temp, μ, 0.5 * (-E1 + ω))
            function integrand3(E1)
                discriminant = (4 * ((E1^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E1^2 * ω^2)
                if discriminant > 0.0
                    return s * pauli_blocking3(E1) / (sqrt(discriminant) * 4π)
                else
                    # @show discriminant, y, E1, s
                    return 0.0
                end
            end
            # integrand3(E1) = real(-(s) * pauli_blocking3(E1) / (sqrt(Complex(4 * ((E1^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E1^2 * ω^2)) * 4π))
            return integrate(integrand3, y, sqrt(param.Λ^2 + m^2))

        elseif s > 4m^2
            y = q * sqrt(1 - 4m^2 / s)
            if y < 1e-5
                return 0.0
            end

            if ω > 0.0
                pauli_blocking1(E2) = numberF(temp, μ, 0.5 * (-E2 - ω)) - numberF(temp, μ, 0.5 * (-E2 + ω))
                function integrand1(E2)
                    discriminant = (4 * ((E2^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E2^2 * ω^2)
                    if discriminant > 0.0
                        return s * pauli_blocking1(E2) / (sqrt(discriminant) * 4π)
                    else
                        # @show discriminant, y, E2, s, s * pauli_blocking1(E2) / sqrt(Complex(discriminant))
                        return 0.0
                    end
                end
                # integrand1(E2) = (s) * pauli_blocking1(E2) / (sqrt((4 * ((E2^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E2^2 * ω^2)) * 4π)
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

function imagpart_phi_q(ω, temp, μ, q, param)
    if q == 0.0
        return imagpart_phi(ω, temp, μ, param)
    end

    m = σ1(temp, μ, param)
    s = ω^2 - q^2
    En_cutoff = sqrt(param.Λ^2 + m^2)
    ω_max = 2*sqrt(param.Λ^2 + q^2/4 + m^2)

    # if s <= -param.Λ^2
    #     return 0.0
    # end

    if ω >= ω_max
        return 0.0
    end

    if (0.0 <= s <= 4m^2)
        return 0.0
    end

    if s < 0.0
        y = q * sqrt(1 - 4m^2 / s)
        if y >= ω_max # what happens when you don't assume it
            return 0.0
        end
        pauli_blocking3(E1) = numberF(temp, μ, 0.5 * (E1 - ω)) - numberF(temp, μ, 0.5 * (E1 + ω)) + numberF(temp, μ, 0.5 * (-E1 - ω)) - numberF(temp, μ, 0.5 * (-E1 + ω))
        function integrand3(E1)
            discriminant = (4 * ((E1^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E1^2 * ω^2)
            if discriminant > 0.0
                return s * pauli_blocking3(E1) / (sqrt(discriminant) * 4π)
            else
                @show discriminant, y, E1, s
                return 0.0
            end
        end
        return -integrate(integrand3, y+1e-5, ω_max)

    else # s > 4m^2
        y = q * sqrt(1 - 4m^2 / s)
        # if y < 1e-5 ## integral is from -y to y
        #     return 0.0
        # end

        if ω > 0.0
            pauli_blocking1(E2) = numberF(temp, μ, 0.5 * (-E2 - ω)) - numberF(temp, μ, 0.5 * (-E2 + ω))
            function integrand1(E2)
                discriminant = ((E2^2 + ω^2 - 4*m^2 - q^2) * q^2 - E2^2 * ω^2)
                if discriminant > 0.0
                    return s * pauli_blocking1(E2) / (sqrt(discriminant) * 4π)
                else
                    @show discriminant, y, E2, s, s * pauli_blocking1(E2) / sqrt(Complex(discriminant))
                    return 0.0
                end
            end
            return integrate(integrand1, -y, y)
        else
            pauli_blocking2(E2) = numberF(temp, μ, 0.5 * (E2 - ω)) - numberF(temp, μ, 0.5 * (E2 + ω))
            integrand2(E2) = (s) * pauli_blocking2(E2) / (sqrt(4 * ((E2^2 + ω^2) / 4 - m^2 - q^2 / 4) * q^2 - E2^2 * ω^2) * 4π)
            return integrate(integrand2, -y, y)
        end
    end
end

function imagpart_phi_q_refactored(ω, temp, μ, q, param)
    if q==0.0
        return imagpart_phi(ω, temp, μ, param)
    end

    # parameters
    m = σ1(temp, μ, param)
    s = ω^2 - q^2
    ω_max = 2*sqrt(param.Λ^2 + q^2/4 + m^2)

    if ω >= ω_max || 0.0<=s<=4m^2
        return 0.0
    end

    if s<0
        return imagpart_phi_q_s_negative(ω, temp, μ, q, m, s, ω_max)
    else # s>4m^2
        return imagpart_phi_q_s_positive(ω, temp, μ, q, m, s, ω_max)
    end
end

function imagpart_phi_q_refactored_k(ω, temp, μ, q, param)
    if q==0.0
        return imagpart_phi(ω, temp, μ, param)
    end

    # parameters
    m = mass_k(temp, μ, q, param)
    s = ω^2 - q^2
    ω_max = 2*sqrt(param.Λ^2 + q^2/4 + m^2)

    if ω >= ω_max || 0.0<=s<=4m^2
        return 0.0
    end

    if s<0
        return imagpart_phi_q_s_negative(ω, temp, μ, q, m, s, ω_max)
    else # s>4m^2
        return imagpart_phi_q_s_positive(ω, temp, μ, q, m, s, ω_max)
    end
end

function imagpart_phi_q_refactored_k(ω, temp, μ, q, m, param)
    if q==0.0
        return imagpart_phi(ω, temp, μ, param)
    end

    # parameters
    s = ω^2 - q^2
    ω_max = 2*sqrt(param.Λ^2 + q^2/4 + m^2)

    if ω >= ω_max || 0.0<=s<=4m^2
        return 0.0
    end

    if s<0
        return imagpart_phi_q_s_negative(ω, temp, μ, q, m, s, ω_max)
    else # s>4m^2
        return imagpart_phi_q_s_positive(ω, temp, μ, q, m, s, ω_max)
    end
end

function imagpart_phi_q_s_negative(ω, temp, μ, q, m, s, ω_max)
    y = q*sqrt(1 - 4m^2/s)
    if y >= ω_max
        return 0.0
    end
    function integrand_1(E1)
        discriminant = (s-4m^2)*q^2 - s*E1^2
        if discriminant > 0.0
            return s*(pauli_blocking_q_1(temp, μ, E1, ω) + pauli_blocking_q_1(temp, μ, -E1, ω)) / (sqrt(discriminant)*4π)
        else
            return 0.0
        end
    end
    return -integrate(integrand_1, y, ω_max)
end

function imagpart_phi_q_s_positive(ω, temp, μ, q, m, s, ω_max)
    y = q*sqrt(1 - 4m^2/s)
    if y < 1e-5
        return 0.0
    end

    function integrand_2(E2)
        discriminant = (s-4m^2)*q^2 - s*E2^2
        if discriminant > 0.0
            return s*pauli_blocking_q_1(temp, μ, E2, ω) / (sqrt(discriminant)*4π)
        else
            return 0.0
        end
    end
    return integrate(integrand_2, -y, y)

end

function pauli_blocking_q_1(temp, μ, x, ω)
    return numberF(temp, μ, 0.5 * (-x - ω)) - numberF(temp, μ, 0.5 * (-x + ω))
end

Π0_phi_q(temp, μ, q, param::Parameters) = Π0_phi(temp, μ, param)

Π0_sigma_q(temp, μ, q, param::Parameters) = Π0_sigma(temp, μ, param)

export imagpart_sigma_q, imagpart_phi_q, Π0_phi_q, Π0_sigma_q
