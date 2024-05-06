"""
    realpart(imagpart, ω, args...)

Returns the real part of the polarisation function given the imaginary part using Krammers-Kronig relations.
"""
function realpart(imagpart::Function, ω, args...)
    integrand(ν) = 2 * ν * imagpart(ν, args...) * (PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2)) / π
    int_sub(ν) = integrand(1 / (1 - ν)) / (1 - ν)^2

    return integrate(integrand, 0.01, 300.0) #+ integrate(int_sub, 0.0, 1.0)
end

function phasesc(imagpart::Function, ω, args...)
    return angle(Complex(realpart(imagpart, ω, args...), -imagpart(ω, args...)))
end

function phaser(imagpart::Function, Π0::Function, ω, args...)
    repi = realpart(imagpart, ω, args...)
    impi = imagpart(ω, args...)
    Π00 = Π0(args...)
    return -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))
end

function phase_tot(imagpart::Function, Π0::Function, ω, args...)
    phase_r = phaser(imagpart, Π0, ω, args...)
    phase_sc = phasesc(imagpart, ω, args...)
    return phase_r + phase_sc
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
