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

function energy_MF(temp,μ,param)
	
end

function number_MF(temp,μ,param)
	
end

export pressure_MF, energy_MF, number_MF