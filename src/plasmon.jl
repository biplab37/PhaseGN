# Polarisation at zero temperature and chemical potential.
function polarisation_imag(ω, q)
    m = σ1(0.01, 0.0, Parameters(κ=0.0))
    s = ω^2 - q^2
    if (s < 0) || (s - 4 * m^2 < 0)
        return 0.0
    else
        p_start = (ω - q * sqrt((s - 4 * m^2) / s)) / 2
        p_end = (ω + q * sqrt((s - 4 * m^2) / s)) / 2
        integrand(Ep) = s / (4 * π * sqrt(4 * (Ep^2 - m^2) * q^2 - (q^2 + ω * (2 * Ep - ω))^2))
        return integrate(integrand, p_start, p_end)
    end
end

# Use kramer-kronig relation to find the real part of the polarisation function from the imaginary part. 
function polarisation_real(ω, q)
    inte(ωp) = 2 * ωp * polarisation_imag(ω, q) / (π * (ωp + ω))
    return PVintegral(inte, 0, 1000, ω, integrate)
end

export polarisation_imag, polarisation_real
