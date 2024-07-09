function Omega(ϕ, temp, μ, param)
    term_vac = -ϕ^2/(2π) + abs(ϕ)^3/(3π)
    term_temp(p) = p*temp*(log(1 + exp(-(sqrt(p^2 + ϕ^2) - μ)/temp)) + log(1 + exp(-(sqrt(p^2 + ϕ^2) + μ)/temp)))/π
    terms = term_vac - integrate(term_temp, 0, param.Λ)
    return terms
end

export Omega
