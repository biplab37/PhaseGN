function realpart(imagpart::Function, ω, args...)
    integrand(ν) = 2 * ν * imagpart(ω, args...) * (PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2)) / π
    int(ν) = integrand(1 / (1 - ν)) / (1 - ν)^2

    return integrate(integrand, 0, 1) + integrate(int, 0, 1)
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

function phase(imagpart::Function, Π0::Function, ω, args...)
    repi = realpart(imagpart, ω, args...)
    impi = imagpart(ω, args...)
    Π00 = Π0(args...)
    phasesc = angle(Complex(repi, -impi))
    phaser = -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))
    return [phasesc, phaser, phasesc + phaser]
end

export realpart, phasesc, phaser, phase