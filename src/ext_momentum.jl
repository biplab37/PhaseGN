# This file contains the code for the case with non zero external momentum

#-------------------------------------------------------------#
#|                      Pseudo Scalar(ϕ)                     |#
#-------------------------------------------------------------#

function get_p1_p2(ω, q, m)

    a = (4m^2 * ω^2 + q^2 * ω^2 - ω^4) / (q^2 - ω^2)
    if a < 0 || a == Inf || a == -Inf
        return 0.0, 0.0
    else
        p1 = abs(q - sqrt(a)) / 2.0
        p2 = (q + sqrt(a)) / 2.0
        return p1, p2
    end
end

"""
	imagpart_ϕ_q(temp,μ,ω,q,param)

Returns the imaginary part of the polarisation function at finite external momentum q for the pseudo scalar channel.
"""
function imagpart_ϕ_q_test(temp, μ, ω, q, param)
    m = σ1(temp, μ, param)
    s = ω^2 - q^2
    E(p) = sqrt(p^2 + m^2)
    p1, p2 = get_p1_p2(ω, q, m)

    integrand1(p) = s * (numberF(temp, μ, -ω + E(p)) - numberF(temp, μ, E(p))) / (8 * π * q * E(p) * sqrt(1 - ((s - 2 * ω * E(p)) / (2 * p * q))^2))
    integrand2(p) = s * (numberF(temp, μ, -ω - E(p)) - numberF(temp, μ, -E(p))) / (8 * π * q * E(p) * sqrt(1 - ((s + 2 * ω * E(p)) / (2 * p * q))^2))
    integrand3(p) = -s * (numberF(temp, μ, -ω + E(p)) - numberF(temp, μ, E(p))) / (8 * π * q * E(p) * sqrt(1 - ((s - 2 * ω * E(p)) / (2 * p * q))^2))
    integrand4(p) = -s * (numberF(temp, μ, -ω - E(p)) - numberF(temp, μ, -E(p))) / (8 * π * q * E(p) * sqrt(1 - ((s + 2 * ω * E(p)) / (2 * p * q))^2))

    ## TODO: - Logic below can be improved
    if ω == 0 || s >= 4(param.Λ^2 + m^2)
        return 0.0
    end
    # s0 = min(2*abs(ω)*m,4m^2)
    # if s>2*abs(ω)*m
    # 	if ω>0
    # 		result = integrate(integrand1,p1,p2)
    # 	else
    # 		result = integrate(integrand2,p1,p2)
    # 	end
    # elseif s<0
    # 	if ω>0
    # 		result = integrate(integrand1,p2,param.Λ)
    # 	else
    # 		result = integrate(integrand2,p2,param.Λ)
    # 	end
    # else
    # 	result = 0.0
    # end
    result = 0.0
    if ω < 0
        if s > 2 * abs(ω) * m
            result += integrate(integrand2, p1, p2) + integrate(integrand4, p1, p2)
        elseif s > 4m^2
            result += integrate(integrand4, p1, p2)
        elseif s < 0
            result += integrate(integrand2, p2, param.Λ) + integrate(integrand3, p1, param.Λ) + integrate(integrand4, p2, param.Λ)
        end
    elseif ω > 0
        if s > 2 * abs(ω) * m
            result += integrate(integrand1, p1, p2) + integrate(integrand3, p1, p2)
        elseif s > 4m^2
            result += integrate(integrand3, p1, p2)
        elseif s < 0
            result += integrate(integrand1, p2, param.Λ) + integrate(integrand3, p2, param.Λ) + integrate(integrand4, p1, param.Λ)
        end
    end

    return result
end
"""
A function which evaluates the delta function numerically and used to test the analytical expression used for further numerical calculation above.
"""
function imagpart_ϕ_q2(temp, μ, ω, q, param)
    m = σ1(temp, μ, param)
    s = ω^2 - q^2
    E(p) = sqrt(p^2 + m^2)
    Ek(p, q, θ) = sqrt(p^2 + q^2 + 2 * p * q * cos(θ) + m^2)

    if ω == 0 || s >= 4(param.Λ^2 + m^2)
        return 0.0
    end
    integ1(p, θ) = s * p * (numberF(temp, μ, -Ek(p, q, θ)) - numberF(temp, μ, E(p))) * DiracDelta(ω - E(p) - Ek(p, q, θ), 0.2) / (8π * E(p) * Ek(p, q, θ))
    integ2(p, θ) = s * p * (numberF(temp, μ, Ek(p, q, θ)) - numberF(temp, μ, -E(p))) * DiracDelta(ω + E(p) + Ek(p, q, θ), 0.2) / (8π * E(p) * Ek(p, q, θ))
    integ3(p, θ) = -s * p * (numberF(temp, μ, Ek(p, q, θ)) - numberF(temp, μ, E(p))) * DiracDelta(ω - E(p) + Ek(p, q, θ), 0.2) / (8π * E(p) * Ek(p, q, θ))
    integ4(p, θ) = -s * p * (numberF(temp, μ, -Ek(p, q, θ)) - numberF(temp, μ, -E(p))) * DiracDelta(ω + E(p) - Ek(p, q, θ), 0.2) / (8π * E(p) * Ek(p, q, θ))

    integrand1(x) = integ1(x[1], x[2])
    integrand2(x) = integ2(x[1], x[2])
    integrand3(x) = integ3(x[1], x[2])
    integrand4(x) = integ4(x[1], x[2])

    if ω > 0
        result = integrate(integrand1, [0.0, 0.0], [param.Λ, 2π]) + integrate(integrand3, [0.0, 0.0], [param.Λ, 2π]) + integrate(integrand4, [0.0, 0.0], [param.Λ, 2π])
    else
        result = integrate(integrand2, [0.0, 0.0], [param.Λ, 2π]) + integrate(integrand3, [0.0, 0.0], [param.Λ, 2π]) + integrate(integrand4, [0.0, 0.0], [param.Λ, 2π])
    end
end

function imagpart_ϕ_q1(temp, μ, ω, q, param)
    β = 1 / temp
    m = σ1(temp, μ, param)
    s = ω^2 - q^2
    Ep(p) = sqrt(p^2 + m^2)

    if s >= 0
        if s <= 4 * m^2 || ω > 2 * param.Λ
            result = 0.0
        else
            p1 = abs(-q + ω * sqrt(1 - 4 * m^2 / s)) / 2
            p2 = (q + ω * sqrt(1 - 4 * m^2 / s)) / 2
            integrand(p) = s * (numberF(temp, μ, -ω + Ep(p)) - numberF(temp, μ, Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s - 2 * ω * Ep(p)) / (2 * p * q))^2))

            result = integrate(integrand, p1, p2)
        end
    else
        p1 = abs(-q + ω * sqrt(1 - 4 * m^2 / s)) / 2
        p2 = (q + ω * sqrt(1 - 4 * m^2 / s)) / 2
        integrand1(p) = -s * (numberF(temp, μ, -ω - Ep(p)) - numberF(temp, μ, -Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s + 2 * ω * Ep(p)) / (2 * p * q))^2))
        result = integrate(integrand1, p1, p2)

        integrand2(p) = -s * (numberF(temp, μ, -ω - Ep(p)) - numberF(temp, μ, -Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s + 2 * ω * Ep(p)) / (2 * p * q))^2)) - s * (numberF(temp, μ, -ω + Ep(p)) - numberF(temp, μ, Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s - 2 * ω * Ep(p)) / (2 * p * q))^2))
        # int1(p) = integrand2(1/(1 - p))/(1 - p)^2

        result += integrate(integrand2, p2, param.Λ)
    end

    return result
end

function imagpart_ϕ_q2(temp, μ, ω, q, m, param)
    β = 1 / temp
    s = ω^2 - q^2
    Ep(p) = sqrt(p^2 + m^2)

    if s >= 0
        if s <= 4 * m^2 || ω > 2 * param.Λ
            result = 0.0
        else
            p1 = abs(-q + ω * sqrt(1 - 4 * m^2 / s)) / 2
            p2 = (q + ω * sqrt(1 - 4 * m^2 / s)) / 2
            integrand(p) = s * (numberF(temp, μ, -ω + Ep(p)) - numberF(temp, μ, Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s - 2 * ω * Ep(p)) / (2 * p * q))^2))

            result = integrate(integrand, p1, p2)
        end
    else
        p1 = abs(-q + ω * sqrt(1 - 4 * m^2 / s)) / 2
        p2 = (q + ω * sqrt(1 - 4 * m^2 / s)) / 2
        integrand1(p) = -s * (numberF(temp, μ, -ω - Ep(p)) - numberF(temp, μ, -Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s + 2 * ω * Ep(p)) / (2 * p * q))^2))
        result = integrate(integrand1, p1, p2)

        integrand2(p) = -s * (numberF(temp, μ, -ω - Ep(p)) - numberF(temp, μ, -Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s + 2 * ω * Ep(p)) / (2 * p * q))^2)) - s * (numberF(temp, μ, -ω + Ep(p)) - numberF(temp, μ, Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s - 2 * ω * Ep(p)) / (2 * p * q))^2))
        int1(p) = integrand2(1 / (1 - p)) / (1 - p)^2

        result += integrate(int1, 1 - 1 / p2, 1)
    end

    return result
end

function imagpart_ϕ_q(temp, μ, ω, q, param)
    β = 1 / temp
    m = σ1(temp, μ, param)
    s = ω^2 - q^2

    integrand1(x) = (numberF(temp, μ, (x - ω) / 2.0) - numberF(temp, μ, (ω + x) / 2.0)) / sqrt(((ω + x)^2 - 4m^2) * q^2 - ω^2 * x^2 - 2 * ω * x * q^2 - q^4)
    integrand2(x) = (numberF(temp, μ, (-x - ω) / 2.0) - numberF(temp, μ, (ω - x) / 2.0)) / sqrt(((-ω + x)^2 - 4m^2) * q^2 - ω^2 * x^2 + 2 * ω * x * q^2 - q^4)
    integrand3(x) = (numberF(temp, μ, (x - ω) / 2.0) - numberF(temp, μ, (ω + x) / 2.0)) / sqrt(((ω + x)^2 - 4m^2) * q^2 - ω^2 * x^2 - 2 * ω * x * q^2 - q^4)
    integrand4(x) = (numberF(temp, μ, (-x - ω) / 2.0) - numberF(temp, μ, (ω - x) / 2.0)) / sqrt(((-ω + x)^2 - 4m^2) * q^2 - ω^2 * x^2 + 2 * ω * x * q^2 - q^4)

    if 0 <= s <= 4m^2
        return 0.0
    else
        y = ω * sqrt(1 - 4 * m^2 / s) / 2.0
        if s < 0.0
            return -s * (integrate(integrand3, y, Inf) + integrate(integrand4, y, Inf)) / 8π
        else
            if ω > 0.0
                return s * (integrate(integrand1, -y, y)) / 8π
            else
                return s * (integrate(integrand2, -y, y)) / 8π
            end
        end
    end
end

"""
	realpart_ϕ_q(temp, μ, ω, q, param)

Returns the real part of the polarisation function at finite external momentum for the pseudo scalar channel.
Uses `imagpart_ϕ_q` and Kramer-Cronig Relation for the calculation
"""
function realpart_ϕ_q(temp, μ, ω, q, param)

    integrand(ν) = 2 * ν * imagpart_ϕ_q(temp, μ, ν, q, param) * (PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2)) / π
    int1(ν) = integrand(1 / (1 - ν)) / (1 - ν)^2

    return integrate(integrand, 0, 1)[1] + integrate(int1, 0, 1)
end

function realpart_ϕ_q(temp, μ, ω, q, m, param)

    integrand(ν) = 2 * ν * imagpart_ϕ_q(temp, μ, ν, q, m, param) * (PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2)) / π
    int1(ν) = integrand(1 / (1 - ν)) / (1 - ν)^2

    return integrate(integrand, 0, 1)[1] + integrate(int1, 0, 1)
end

@doc raw"""
	phasesc_ϕ_q(temp,μ,ω,param)

Returns the scattering part of the phase for ϕ₁ and ϕ₂ at external momentum q.

Also consider the equivalent function 

	phasesc_ϕ_q(temp,μ,ω,m,param)

where you supply the values of $\bar{\sigma}_1 = m$. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phasesc_ϕ_q(temp, μ, ω, q, param)
    repi = realpart_ϕ_q(temp, μ, ω, q, param)
    impi = imagpart_ϕ_q(temp, μ, ω, q, param)
    return angle(Complex(repi, -impi))
end
function phasesc_ϕ_q(temp, μ, ω, q, m, param)
    repi = realpart_ϕ_q(temp, μ, ω, q, m, param)
    impi = imagpart_ϕ_q(temp, μ, ω, q, m, param)

    return angle(Complex(repi, -impi))
end

@doc raw"""
	phaser_ϕ_q(temp,μ,ω,param)

Returns the resonant part of the phase for ϕ₁ and ϕ₂ at external momentum q.

Also consider the equivalent function 

	phaser_ϕ_q(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency 
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases
at different vaules of frequencies at a fixed temp and μ.
"""
function phaser_ϕ_q(temp, μ, ω, q, param)
    repi = realpart_ϕ_q(temp, μ, ω, q, param)
    impi = imagpart_ϕ_q(temp, μ, ω, q, param)
    Π00 = Π0_ϕ(temp, μ, param)

    return -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))
end

function phaser_ϕ_q(temp, μ, ω, q, m, Π00, param)
    repi = realpart_ϕ_q(temp, μ, ω, q, m, param)
    impi = imagpart_ϕ_q(temp, μ, ω, q, m, param)

    return -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))
end

@doc raw"""
	phase_ϕ_q(temp,μ,ω,param)

Returns all the phases for ϕ₁ and ϕ₂ at external momentum q in an array 
	[scattered phase, resonant phase, total phase]

Also consider the equivalent function 

	phase_ϕ_q(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phase_ϕ_q(temp, μ, ω, q, param)
    repi = realpart_ϕ_q(temp, μ, ω, q, param)
    impi = imagpart_ϕ_q(temp, μ, ω, q, param)
    Π00 = Π0_ϕ(temp, μ, param)

    phasesc = angle(Complex(repi, -impi))
    phaser = -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))

    return [phasesc, phaser, phasesc + phaser]
end

function phase_ϕ_q(temp, μ, ω, q, m, Π00, param)
    repi = realpart_ϕ_q(temp, μ, ω, q, m, param)
    impi = imagpart_ϕ_q(temp, μ, ω, q, m, param)

    phasesc = angle(Complex(repi, -impi))
    phaser = -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))

    return [phasesc, phaser, phasesc + phaser]
end

#          ╭──────────────────────────────────────────────────────────╮
#          │                        Scalar(σ)                         │
#          ╰──────────────────────────────────────────────────────────╯

"""
	imagpart_σ_q(temp, μ, ω,param)

Returns the imaginary part of the polarisation for σ₁ and σ₂ at external momentum q.
"""
function imagpart_σ_q(temp, μ, ω, q, param)
    β = 1 / temp
    m = σ1(temp, μ, param)
    s = ω^2 - q^2
    Ep(p) = sqrt(p^2 + m^2)

    if s >= 0
        if s <= 4 * m^2 || ω > 2 * param.Λ
            result = 0.0
        else
            p1 = abs(-q + ω * sqrt(1 - 4 * m^2 / s)) / 2
            p2 = (q + ω * sqrt(1 - 4 * m^2 / s)) / 2
            integrand(p) = (s - 4 * m^2) * (numberF(temp, μ, -ω + Ep(p)) - numberF(temp, μ, Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s - 2 * ω * Ep(p)) / (2 * p * q))^2))

            result = integrate(integrand, p1, p2)
        end
    else
        p1 = abs(-q + ω * sqrt(1 - 4 * m^2 / s)) / 2
        p2 = (q + ω * sqrt(1 - 4 * m^2 / s)) / 2
        integrand1(p) = -(s - 4 * m^2) * (numberF(temp, μ, -ω - Ep(p)) - numberF(temp, μ, -Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s + 2 * ω * Ep(p)) / (2 * p * q))^2))
        result = integrate(integrand1, p1, p2)

        integrand2(p) = -(s - 4 * m^2) * (numberF(temp, μ, -ω - Ep(p)) - numberF(temp, μ, -Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s + 2 * ω * Ep(p)) / (2 * p * q))^2)) - (s - 4 * m^2) * (numberF(temp, μ, -ω + Ep(p)) - numberF(temp, μ, Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s - 2 * ω * Ep(p)) / (2 * p * q))^2))
        result += integrate(integrand2, p2, param.Λ)
    end

    return result
end

function imagpart_σ_q(temp, μ, ω, q, m, param)
    β = 1 / temp
    s = ω^2 - q^2
    Ep(p) = sqrt(p^2 + m^2)

    if s >= 0
        if s <= 4 * m^2 || ω > 2 * param.Λ
            result = 0.0
        else
            p1 = abs(-q + ω * sqrt(1 - 4 * m^2 / s)) / 2
            p2 = (q + ω * sqrt(1 - 4 * m^2 / s)) / 2
            integrand(p) = (s - 4 * m^2) * (numberF(temp, μ, -ω + Ep(p)) - numberF(temp, μ, Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s - 2 * ω * Ep(p)) / (2 * p * q))^2))

            result = integrate(integrand, p1, p2)
        end
    else
        p1 = abs(-q + ω * sqrt(1 - 4 * m^2 / s)) / 2
        p2 = (q + ω * sqrt(1 - 4 * m^2 / s)) / 2
        integrand1(p) = -(s - 4 * m^2) * (numberF(temp, μ, -ω - Ep(p)) - numberF(temp, μ, -Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s + 2 * ω * Ep(p)) / (2 * p * q))^2))
        result = integrate(integrand1, p1, p2)

        integrand2(p) = -(s - 4 * m^2) * (numberF(temp, μ, -ω - Ep(p)) - numberF(temp, μ, -Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s + 2 * ω * Ep(p)) / (2 * p * q))^2)) - (s - 4 * m^2) * (numberF(temp, μ, -ω + Ep(p)) - numberF(temp, μ, Ep(p))) / (8 * π * q * Ep(p) * sqrt(1 - ((s - 2 * ω * Ep(p)) / (2 * p * q))^2))
        result += integrate(integrand2, p2, param.Λ)
    end

    return result
end

"""
	realpart_σ_q(temp, μ, ω,param)

Returns the real part of the polarisation for σ₁ and σ₂ at external momentum q.
"""
function realpart_σ_q(temp, μ, ω, q, param)

    integrand(ν) = 2 * ν * imagpart_σ_q(temp, μ, ν, q, param) * (PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2)) / π
    int(ν) = integrand(1 / (1 - ν)) / (1 - ν)^2

    return integrate(integrand, 0, 1, maxevals=10000)[1] + integrate(int, 0, 1)
end

function realpart_σ_q(temp, μ, ω, q, m, param)
    integrand(ν) = 2 * ν * imagpart_σ_q(temp, μ, ν, q, m, param) * (PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2))
    int(ν) = integrand(1 / (1 - ν)) / (1 - ν)^2

    return integrate(integrand, 0, 1, maxevals=10000)[1] + integrate(int, 0, 1)
end

@doc raw"""
	phasesc_σ_q(temp,μ,ω,param)

Returns the scattering part of the phase for σ₁ and σ₂ at external momentum q.

Also consider the equivalent function 

    phasesc_ϕ_q(temp,μ,ω,m,param)

where you supply the values of $\bar{\sigma}_1 = m$. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phasesc_σ_q(temp, μ, ω, q, param)
    repi = realpart_σ_q(temp, μ, ω, q, param)
    impi = imagpart_σ_q(temp, μ, ω, q, param)
    return angle(Complex(repi, -impi))
end

function phasesc_σ_q(temp, μ, ω, q, m, param)
    repi = realpart_σ_q(temp, μ, ω, q, m, param)
    impi = imagpart_σ_q(temp, μ, ω, q, m, param)

    return angle(Complex(repi, -impi))
end

@doc raw"""
	phaser_σ_q(temp,μ,ω,param)

Returns the resonant part of the phase for σ₁ and σ₂ at external momentum q.

Also consider the equivalent function 

    phaser_ϕ_q(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency 
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases
at different vaules of frequencies at a fixed temp and μ.
"""
function phaser_σ_q(temp, μ, ω, q, param)
    repi = realpart_σ_q(temp, μ, ω, q, param)
    impi = imagpart_σ_q(temp, μ, ω, q, param)
    Π00 = Π0_σ(temp, μ, param)

    return -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))
end

function phaser_σ_q(temp, μ, ω, q, m, Π00, param)
    repi = realpart_σ_q(temp, μ, ω, q, m, param)
    impi = imagpart_σ_q(temp, μ, ω, q, m, param)

    return -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))
end

@doc raw"""
	phase_σ_q(temp,μ,ω,param)

Returns all the phases for σ₁ and σ₂ at external momentum q in an array 

	[scattered phase, resonant phase, total phase]

Also consider the equivalent function 

    phase_σ_q(temp,μ,ω,m,Π00,param)

where you supply the values of $\bar{\sigma}_1 = m$ and $\Pi 00$ the frequency
and momentum independent part of the polarisation. Since these values only depend 
on the temp and μ, this function is more efficient if you want to calculate phases 
at different vaules of frequencies at a fixed temp and μ.
"""
function phase_σ_q(temp, μ, ω, q, param)
    repi = realpart_σ_q(temp, μ, ω, q, param)
    impi = imagpart_σ_q(temp, μ, ω, q, param)
    Π00 = Π0_σ(temp, μ, param)

    phasesc = angle(Complex(repi, -impi))
    phaser = -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))

    return [phasesc, phaser, phasesc + phaser]
end

function phase_σ_q(temp, μ, ω, q, m, Π00, param)
    repi = realpart_σ_q(temp, μ, ω, q, m, param)
    impi = imagpart_σ_q(temp, μ, ω, q, m, param)

    phasesc = angle(Complex(repi, -impi))
    phaser = -angle(Complex(repi - (repi^2 + impi^2) / Π00, -impi))

    return [phasesc, phaser, phasesc + phaser]
end

export imagpart_ϕ_q, imagpart_σ_q, realpart_ϕ_q, realpart_σ_q, phasesc_ϕ_q, phasesc_σ_q, phaser_ϕ_q, phaser_σ_q, phase_ϕ_q, phase_σ_q
