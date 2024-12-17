"""
    realpart(imagpart, ω, args...)

Returns the real part of the polarisation function given the imaginary part using Krammers-Kronig relations.
"""
function realpart(imagpart::Function, ω, T, μ, param)
    integrand(ν) = 2 * ν * imagpart(ν, T, μ, param) * (PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2)) / π
    int_sub(ν) = integrand(1 / (1 - ν)) / (1 - ν)^2

    return integrate(integrand, 0.001, param.Λ)# + integrate(int_sub, 0.0, 1.0)
end

function realpart(imagpart::Function, ω, T, μ, q, param)
    m = σ1(T, μ, param)
    cutoff = 2*sqrt(param.Λ^2 + q^2/4 + m^2)
    integrand(ν) = 2 * ν * (imagpart(ν, T, μ, q, m, param) * PrincipalValue(ν^2 - ω^2) - imagpart(ν, T, μ, 0.0, param) * PrincipalValue(ν^2)) / π
    return integrate(integrand, 0.0, cutoff)
end

function realpart(imagpart::Function, ω, T, μ, q, m, param::Parameters)
    cutoff = 2*sqrt(param.Λ^2 + q^2/4 + m^2)
    integrand(ν) = 2 * ν * (imagpart(ν, T, μ, q, m, param) * PrincipalValue(ν^2 - ω^2) - imagpart(ν, T, μ, 0.0, m, param) * PrincipalValue(ν^2)) / π
    return integrate(integrand, 0.0, cutoff)
end

function realpart_k(imagpart::Function, ω, T, μ, q, param)
    m = mass_k(T, μ, q, param)
    cutoff = 2*sqrt(param.Λ^2 + q^2/4 + m^2)
    integrand(ν) = 2 * ν * (imagpart(ν, T, μ, q, m, param) * PrincipalValue(ν^2 - ω^2) - imagpart(ν, T, μ, 0.0, param) * PrincipalValue(ν^2)) / π
    return integrate(integrand, 0.0, cutoff)
end

function realpart_3_k(imagpart::Function, ω, T, μ, q, param)
    m = mass_k(T, μ, q, param)
    integrand(ν) = 2 * ν * imagpart(ν, T, μ, q, m, param) / (π * (ν + ω))
    return (param.Λ - 1) / π - PVintegral(integrand, 0.0, 2.1 * param.Λ, ω, integrate)
end

function realpart_3(imagpart::Function, ω, T, μ, q, m, param)
    integrand(ν) = 2 * ν * imagpart(ν, T, μ, q, m, param) / (π * (ν + ω))
    return (param.Λ - 1) / π - PVintegral(integrand, 0.0, sqrt(5)* param.Λ, ω, integrate)
end

function realpart_2(imagpart::Function, ω, T, μ, q, param)
    integrand(ν) = 2 * ν * imagpart(ν, T, μ, q, param) * PrincipalValue(ν^2 - ω^2) / π
    return (param.Λ - 1) / π - integrate(integrand, 0.0, 2.1 * param.Λ)
end

"""
    phasesc(imagpart::Function, ω, args...)

"""
function phasesc(imagpart::Function, ω, args...)
    return angle(Complex(realpart(imagpart, ω, args...), -imagpart(ω, args...)))
end

function phaser(imagpart::Function, Π0::Function, ω, args...)
    repi = realpart(imagpart, ω, args...)
    impi = imagpart(ω, args...)
    Π00 = Π0(args...)
    return -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))
end

"""
    phase_tot(imagpart::Function, Π0::Function, ω, args...)

Returns the total phase for a given imaginary part of the function.
"""
function phase_tot(imagpart::Function, Π0::Function, ω, args...)
    phase_r = phaser(imagpart, Π0, ω, args...)
    phase_sc = phasesc(imagpart, ω, args...)
    return phase_r + phase_sc
end

function phasetot(imagpart::Function, ω, args...)
    repi = realpart_3(imagpart, ω, args...)
    impi = imagpart(ω, args...)
    return angle(Complex(repi, impi))
end
"""
    phase(imagpart, Π0, ω, args...)

Returns the phases for a given frequency ω, temperature temp, chemical potential μ, and momentum q. One needs to provide the imaginary part of the polarisation function imagpart, the momenta independent real part of the polarisation function Π0.
The returned values are respectively the scattering phase, resonant phase and the total phase.
If one inputs an array of frequencies ωrange, the function will return an array of the phases for each frequency.
"""
function phase(imagpart::Function, Π0::Function, ω, args...)
    repi = realpart(imagpart, ω, args...)
    impi = imagpart(ω, args...)
    Π00 = Π0(args...)
    phasesc = angle(Complex(repi, -impi))
    phaser = -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))
    return [phasesc, phaser, phasesc + phaser]
end

function phase(imagpart::Function, Π0::Function, ωrange::AbstractArray, args...)
    phases = zeros(length(ωrange), 3)
    for (i, ω) in enumerate(ωrange)
        phases[i, :] .= phase(imagpart, Π0, ω, args...)
    end
    return phases
end

export realpart, phasesc, phaser, phase_tot, phase
