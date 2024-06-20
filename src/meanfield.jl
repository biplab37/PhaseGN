# This file contains the code to calculate quantities at Mean Field approximation.

"""
    pressure_MF(temp, μ, param::Parameters; norm)

Returns the pressure in the mean field approximation at a given temperature and chemical potential μ.
"""
function pressure_MF(temp, μ, param::Parameters; norm=false::Bool)
    σ₁ = σ1(temp, μ, param)
    M = param.M
    temp_independent = (M * σ₁^2 / 2π) - (σ₁^3 / 3π) - (M^3 / 6π)
    temp_dependent = -temp^3 * (σ₁ * (reli2(-exp(-(σ₁ - μ) / temp)) + reli2(-exp(-(σ₁ + μ) / temp))) / temp + reli3(-exp(-(σ₁ - μ) / temp)) + reli3(-exp(-(σ₁ + μ) / temp))) / π
    if norm
        return (temp_independent + temp_dependent) / temp^3
    else
        return temp_independent + temp_dependent
    end
end

#TODO: finish the implementation of the mean field approximation at zero chemical potential
function pressure_MF(temp, param::Parameters; norm=false::Bool)
    return pressure_MF(temp, 0.0, param, norm=norm)
end

function pressure_MF(trange::AbstractVector, μ, param::Parameters, norm=false)
    pres = zeros(length(trange))
    for (i,T) in enumerate(trange)
        pres[i] = pressure_MF(T, μ, param, norm=norm)
    end
    return pres
end

"""
    energy_MF(temp, μ, param::Parameters; norm)

Returns the energy in the mean field approximation at a given temperature and chemical potential μ.
"""
function energy_MF(temp, μ, param::Parameters; norm=false::Bool)
    β = 1 / temp
    σ = σ1(temp, μ, param)
    M = σ1(0.01, μ, param)
    energy = β^3 * (-3 * M * σ^2 + 2 * σ^3 + M^3) / (6 * pi) + (β * σ * (β * σ * (log(1 + exp(-β * (σ - μ))) + log(1 + exp(-β * (σ + μ)))) - 2 * reli2(-exp(-β * (σ - μ))) - 2 * reli2(-exp(-β * (σ + μ)))) - 2 * (reli3(-exp(-β * (σ - μ))) + reli3(-exp(-β * (σ + μ))))) / pi

    if norm
        return energy
    else
        return energy * temp^3
    end
end

function energy_MF(trange::AbstractVector, μ, param::Parameters, norm=false::Bool)
    ener = zeros(length(trange))
    for (i, T) in enumerate(trange)
        ener[i] = energy_MF(T, μ, param, norm = norm)
    end
end

"""
    number_MF(temp, μ, param::Parameters; norm)

Returns the number in the mean field approximation at a given temperature and chemical potential μ.
"""
function number_MF(temp, μ, param::Parameters; norm=false::Bool)
    σ = σ1(temp, μ, param)
    β = 1 / temp
    number = temp * (β * σ * (log(1 + exp(-β * (σ - μ))) - log(1 + exp(-β * (σ + μ)))) - reli2(-exp(-β * (σ - μ))) + reli2(-exp(-β * (σ + μ)))) / (π)

    if norm
        return number / temp^2
    else
        return number
    end
end

function number_MF(trange::AbstractVector, μ, param::Parameters, norm=false::Bool)
    numb = zeros(length(trange))
    for (i, T) in enumerate(trange)
        numb[i] = number_MF(T, μ, param, norm = norm)
    end
end

export pressure_MF, energy_MF, number_MF
