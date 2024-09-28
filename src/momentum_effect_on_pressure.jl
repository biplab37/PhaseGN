## This file contains code to quantify the effect of the external momentum to the pressure.

"""
    _integrand_boosted(temp, mu, ω, q, param::Parameters)

This function returns the integrand for calculating the pressure with boosted value of the
stationary phase shift as a function of frequency and external momentum.
"""
function _integrand_boosted_phi(temp, mu, ω, q, param::Parameters)
    factor = 1.0 #calculate and replace this value
    phase_shift = phase_phi(ω, temp, mu, param)[3]
    return factor * phase_shift
end

"""
    _integrand_momentum_dep(temp, mu, ω, q, param::Parameters)

This function returns the integrand for calculating the pressure with boosted value of the
stationary phase shift as a function of frequency and external momentum.
"""
function _integrand_momentum_dep(temp, mu, ω, q, param::Parameters)
    factor = 1.0 #calculate and replace this value
    phase_shift = phasetot_phi(temp, mu, ω, q, param)[3]
    return factor * phase_shift
end
