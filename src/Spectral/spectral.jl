
## The following function and the related helper functions calculates the imaginary
## part of the spectral function from the imaginary part of the dielectric function
## using the formula presented in 10.1103/PhysRevB.77.081411.

function Enkq(T, mu, k, q, θ, s, param::Parameters)
    # s = band index ±1
    new_k2 = k^2 + q^2 + 2 * k * q * cos(θ)
    m = mass_k(T, mu, sqrt(new_k2), param)
    return m, (s * sqrt(new_k2 + m^2))
    # Fermi surface at k=0, s positive and negative for conduction and valence band.
end

function polini_integrand(T, mu, ω, k, q, θ, s, sp, param::Parameters)
    m, Ensp = Enkq(T, mu, k, q, θ, sp, param)
    sign = 1
    if ω > 0
        if !(0 < Ensp < ω)
            return 0.0
        end
        sign = 1
    else
        if !(ω < Ensp < 0.0)
            return 0.0
        end
        sign = -1
    end
    #TODO: this should not be imaginary part of the polarisation function but the inverse of it.
    invdielec(Q, omega) = imagpart_phi_q_refactored_m(omega, T, mu, Q, m, param)

    return q*(2*π/(param.Λ - 1.0))*invdielec(q, ω - Ensp)*(1 + s*sp*cos(θ))*sign/(8*π^2)
end

function spectral_polini(T, mu, ω, k, s, param::Parameters)
    f(x, sp) = polini_integrand(T, mu, ω, k, x[1], x[2], s, sp, param)
    intt(f::Function) = integrate(f, [0.0, 0.0], [param.Λ, 2π])
    return intt(x -> f(x, +1)) + intt(x -> f(x, -1))
end