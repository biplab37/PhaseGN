"""
    bisection(func::Function,start::Number,finish::Number,iteration::Integer=20)
Finds the zero of a function `func` inside a given interval (`start`,`finish`). Might fail 
if there are multiple zeros or no zeros inside the interval.
"""
function bisection(func::Function, start::Number, finish::Number, iteration::Integer=20)
	mid = (start + finish)/2.0

    for i in 1:iteration
        if func(mid)*func(finish) > 0
            finish = mid
        else
        	start = mid
        end
    	mid = (start + finish)/2.0
    end
    return mid
end

"""
	DiracDelta(input::Float64,δ::Float64 = 1e-3)
Dirac Delta Function that can be used in a numerical integration.
"""

function DiracDelta(input, δ = 0.05)
	return exp(-input^2/(4*δ))/(2*sqrt(pi*δ))
end

"""
	PrincipalValue(x::Float64,ϵ::Float64 = 1e-3)
Useful function in Principle Value Integration.
"""
function PrincipalValue(x,ϵ = 1e-3)
	if abs(x)<ϵ
		return 0.0
	else
		return 1/x
	end
end

"""
	numberF(temp::Float64, μ::Float64, Energy::Float64)
Returns the number density at a given `Energy`, chemical potential `μ` and `temp` for Fermionic Species.
"""
function numberF(temp, μ, Energy)
	β = 1/temp
    return 1/(1+exp(β*(Energy - μ)))
end

# using QuadGK
function integrate(func::Function,start::Number, finish::Number;maxevals=10000)
	return hquadrature(func,start,finish,reltol=1e-3,maxevals=maxevals)[1]
	# return quadgk(func,start,finish,rtol=1e-3,maxevals=maxevals)[1]
end

function integrate(func,start::Vector, finish::Vector;maxevals=100000)
	return hcubature(func,start,finish,reltol=1e-3,maxevals=maxevals)[1]
end