"""
	pressure(func,temp,μ,param::Parameters)

Returns the pressure corrosponding to the phase function `func(temp,μ,ω,param)` 
"""
function pressure(func::Function,temp,μ,param::Parameters)
	integrand(s) = -(log(exp(sqrt(s)/temp)-1)-sqrt(s)/temp)/(2*π^2*temp^2)*func(temp,μ,sqrt(s),param)
	int1(s) = integrand(1/(1-s))/(1-s)^2
	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]+ hquadrature(int1,0,1,abstol=1e-2,maxevals=10000)[1]
end

function pressure(func::Function,temp,μ,m,param::Parameters)
	integrand(s) = -(log(exp(sqrt(s)/temp)-1)-sqrt(s)/temp)/(2*π^2*temp^2)*func(temp,μ,sqrt(s),m,param)
	int1(s) = integrand(1/(1-s))/(1-s)^2
	return hquadrature(integrand,0,1,reltol=1e-2,maxevals=10000)[1]+ hquadrature(int1,0,1,abstol=1e-2,maxevals=10000)[1]
end

function pressure(func::Function,temp,μ,m,M,Γ,param::Parameters)
	integrand(s) = -(log(exp(sqrt(s)/temp)-1)-sqrt(s)/temp)/(2*π^2*temp^2)*func(temp,μ,sqrt(s)m,M,Γ,param)
	int1(s) = integrand(1/(1-s))/(1-s)^2
	return hquadrature(integrand,0,1,reltol=1e-3,maxevals=10000)[1]+ hquadrature(int1,0,1,abstol=1e-3,maxevals=10000)[1]
end

export pressure
