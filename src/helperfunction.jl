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
function integrate(func::Function,start, finish;maxevals=10000)
	return hquadrature(func,start,finish,reltol=1e-3,maxevals=maxevals)[1]
	# return quadgk(func,start,finish,rtol=1e-3,maxevals=maxevals)[1]
end