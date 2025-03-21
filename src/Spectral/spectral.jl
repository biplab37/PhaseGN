
## The following function and the related helper functions calculates the imaginary
## part of the spectral function from the imaginary part of the dielectric function
## using the formula presented in 10.1103/PhysRevB.77.081411.

function Enkq(k, q, θ, s)
    # s = band index ±1
    new_k2 = k^2 + q^2 + 2 * k * q * cos(θ)
    m = mass_k(T, mu, sqrt(new_k2), param)
    return s * sqrt(new_k2 + m^2) - m
end

function polini_integrand(T, mu, ω, k, q, θ, s, sp, param::Parameters)
    Ensp = Enkq(k, q, θ, sp)
    if ω > 0
        if !(0 < Ensp < ω)
            return 0.0
        end
    else
        if !(ω < Ensp < 0.0)
            return 0.0
        end
    end
end

function spectral_polini(T, mu, ω, k, s, param::Parameters)
    f(x, sp) = polini_integrand(T, mu, ω, k, x[1], x[2], s, sp, param)
    intt(f::Function) = integrate2D(f, [0.0, 0.0], [param.Λ, 2π])
    return intt(x -> f(x, +1)) + intt(x -> f(x, -1))
end