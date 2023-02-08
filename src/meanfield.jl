# This file contains the code to calculate quantities at Mean Field approximation

function pressure_MF(temp,μ,param::Parameters;norm=false::Bool)
	σ₁ = σ1(temp,μ,param)
    temp_independent = param.M*σ₁^2/2π - σ₁^3/3π - param.M^3/6π
    temp_dependent = -temp^3*(σ₁*(reli2(-exp(-(σ₁ - μ)/temp)) + reli2(-exp(-(σ₁ + μ)/temp)))/temp + reli3(-exp(-(σ₁ - μ)/temp)) + reli3(-exp(-(σ₁ + μ)/temp)))/π
    
    if norm
        return (temp_independent + temp_dependent)/temp^3
    else
        return temp_independent + temp_dependent
    end
end

function energy_MF(temp,μ,param::Parameters;norm=false::Bool)
	β = 1/temp
	σ = σ1(temp,μ,param)
	energy = β^3*(-3*param.M*σ^2 + 2*σ^3 + param.M^3)/(6*pi) + (β*σ*(β*σ*(log(1 + exp(-β*(σ - μ))) + log(1 + exp(-β*(σ + μ)))) - 2*reli2(-exp(-β*(σ - μ))) - 2*reli2(-exp(-β*(σ + μ)))) - 2*(reli3(-exp(-β*(σ - μ))) + reli3(-exp(-β*(σ + μ)))))/pi
    
    if norm
        return energy
    else
        return energy*temp^3
    end
end

function number_MF(temp,μ,param::Parameters;norm=false::Bool)
    σ = σ1(temp,μ,param)
    β = 1/temp
    number = temp*(β*σ*(log(1 + exp(-β*(σ - μ))) - log(1 + exp(-β*(σ + μ)))) - reli2(-exp(-β*(σ - μ)) ) + reli2(-exp(-β*(σ + μ)) ) )/(π)

    if norm
        return number/temp^2
    else
        return number
    end
end

export pressure_MF, energy_MF, number_MF