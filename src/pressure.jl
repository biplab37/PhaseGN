"""
	pressure(func,temp,μ,κ)

Returns the pressure corrosponding to the phase function `func(temp,μ,ω,κ)` 
"""
function pressure(func,temp,μ,κ)
	integrand(s) = -(log(exp(sqrt(s)/temp)-1)-sqrt(s)/temp)/(2*π^2*temp^2)*func(temp,μ,sqrt(s),κ)
	int1(s) = integrand(1/(1-s))/(1-s)^2
	return hquadrature(integrand,0,1,reltol=1e-3,maxevals=10000)[1]+ hquadrature(int1,0,1,abstol=1e-3,maxevals=10000)[1]
end

export pressure
