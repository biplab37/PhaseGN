## This file calculates the dynamic structure factor of the system given by the formula
## Sₘ(ω,q) = 1/(2π) (-2Im(Dₘ(ω,q)))/(exp(ω/T)-1)
## where spectral function is given by ρₘ(ω,q) = - 2Im(Dₘ(ω,q))

function spectral_function_ϕ(T, μ, ω, q, param::Parameters)
    return -2 * imagpart_ϕ_q(T, μ, ω, q, param)
end

function spectral_function_σ(T, μ, ω, q, param::Parameters)
    return -2 * imagpart_σ_q(T, μ, ω, q, param)
end

function structure_factor_ϕ(T, μ, ω, q, param::Parameters)
    return 1 / (2π) * (-2 * imagpart_ϕ_q(T, μ, ω, q, param)) / (exp(ω / T) - 1)
end

function structure_factor_σ(T, μ, ω, q, param::Parameters)
    return 1 / (2π) * (-2 * imagpart_σ_q(T, μ, ω, q, param)) / (exp(ω / T) - 1)
end

export structure_factor_ϕ, structure_factor_σ
