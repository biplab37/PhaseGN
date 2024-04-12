function realpart(imagpart::Function, ω, args...)
    integrand(ν) = 2 * ν * imagpart(ν, args...) * (PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2)) / π
    int_sub(ν) = integrand(1 / (1 - ν)) / (1 - ν)^2

    return integrate(integrand, 0, 1) + integrate(int_sub, 0, 1)
end

function phasesc(imagpart::Function, ω, args...)
    return angle(Complex(realpart(imagpart, ω, args...), imagpart(ω, args...)))
end

function phaser(imagpart::Function, Π0::Function, ω, args...)
    repi = realpart(imagpart, ω, args...)
    impi = imagpart(ω, args...)
    Π00 = Π0(args...)
    return -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))
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
    phaser = zero(length(ωrange))
    phasesc = zero(length(ωrange))
    phase = zero(length(ωrange))
    for (i, ω) in enumerate(ωrange)
        phaser[i], phasesc[i], phase[i] = phase(imagpart, Π0, ω, args...)
    end
    return [phasesc, phaser, phase]
end

export realpart, phasesc, phaser, phase
