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
	PrincipleValue(x::Float64,ϵ::Float64 = 1e-3)
Useful function in Principle Value Integration.
"""
function PrincipleValue(x,ϵ = 1e-3)
	if abs(x)<ϵ
		return 0.0
	else
		return 1/x
	end
end