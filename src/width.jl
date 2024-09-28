# This file contains code for the calculation of width

function realpartI2(temp, μ, param)
    m = σ1(temp, μ, param)
    M = M_phi(temp, μ, param)
    Ep(p) = sqrt(p^2 + m^2)
    integrand(p) = p * (1 - numberF(temp, -μ, Ep(p)) - numberF(temp, μ, Ep(p))) * (PrincipalValue(2 * Ep(p) - M, 1e-3) + 1 / (M + 2 * Ep(p))) / (8 * π * Ep(p)^2)
    return integrate(integrand, 0.0, param.Λ)
end

function imagpartI2(temp, μ, param)
    m = σ1(temp, μ, param)
    M = M_phi(temp, μ, param)
    Nσ = M^2 - 4 * m^2
    if Nσ > 0 && abs(M) < 2 * sqrt(param.Λ^2 + m^2)
        return (1 - numberF(temp, -μ, M / 2) - numberF(temp, μ, M / 2)) / (8 * M)
    else
        return 0.0
    end
end

function I1(temp, μ, param)
    m = σ1(temp, μ, param)
    Ep(p) = sqrt(p^2 + m^2)
    integrand(p) = 1 / (2 * π) - 2 * p * (1 - numberF(temp, -μ, Ep(p)) - numberF(temp, μ, Ep(p))) / (4 * π * Ep(p))
    return -1 / (2 * π) + integrate(integrand, 0.0, param.Λ)
end

function Γϕ(temp, μ, param)
    M = M_phi(temp, μ, param)
    imp = imagpartI2(temp, μ, param)
    rep = realpartI2(temp, μ, param)
    i1 = I1(temp, μ, param)

    Γ = sqrt(i1 * imp / ((imp^2 + rep^2) * M))
    return Γ
end

function Γϕ2(temp, μ, param)
    m = σ1(temp, μ, param)
    M = M_phi(temp, μ, param)
    imp = imagpartI2(temp, μ, param)
    rep = realpartI2(temp, μ, param)
    i1 = I1(temp, μ, param)
    x = M^2 - (i1 * rep / ((imp^2 + rep^2)))

    if x >= 0
        return 2 * sqrt(x)
    else
        return 0.0
    end
end

function phaseBW_ϕ(temp, μ, ω, param)
    Γ = Γϕ(temp, μ, param)
    m = σ1(temp, μ, param)
    M = M_phi(temp, μ, param)
    if 2 * m >= M
        if ω <= 2 * m
            return 0.0
        else
            return π
        end
    else
        return π * (atan((ω^2 - M^2) / (M * Γ)) - atan((4 * m^2 - M^2) / (M * Γ))) / (π / 2 - atan((4 * m^2 - M^2) / (M * Γ)))
    end
end

function phaseBW_ϕ(temp, μ, ω, m, M, Γ, param)
    if Γ == 0
        if ω < 2 * m
            return 0.0
        else
            return π
        end
    else
        return π * (atan((ω^2 - M^2) / (M * Γ)) - atan((4 * m^2 - M^2) / (M * Γ))) / (π / 2 - atan((4 * m^2 - M^2) / (M * Γ)))
    end
end

function Γ_ϕ(temp, μ, param)
    dω = 1e-2

    ωX = M_phi(temp, μ, param)
    m = σ1(temp, μ, param)

    impi = imagpart_phi(temp, μ, ωX, m, param)
    rep(ω) = fullrealpart_phi(temp, μ, ω, m, param) - realpart_phi(temp, μ, ω, m, param)
    repi = rep(ωX)

    drepi = (rep(sqrt(ωX^2 + dω)) - rep(ωX)) / dω
    dimpi = 0.0 #( imagpart_ϕ(temp,μ,ωX+dω,m,param) - impi )/dω


    return 2 * impi / ((drepi - 2 * repi * (repi * drepi + impi * dimpi) / (repi^2 + impi^2)))
end

export Γϕ, Γϕ2, phaseBW_ϕ, Γ_ϕ
