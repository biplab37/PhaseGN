## Contains code about the calculation of pressure using the refactored version used for the paper.

function delta_integrand(func::Function, ω, q, T, mu, param::Parameters)
    return q * func(ω, q, T, mu, param) / ((exp(ω / T) - 1) * 2π^2)
end
