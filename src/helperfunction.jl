# using QuadGK
function integrate(func::Function,start::Number, finish::Number;maxevals=10000)
	return hquadrature(func,start,finish,reltol=1e-3,maxevals=maxevals)[1]
	# return quadgk(func,start,finish,rtol=1e-3,maxevals=maxevals)[1]
end

function integrate(func,start::Vector, finish::Vector;maxevals=100000)
	return hcubature(func,start,finish,reltol=1e-3,maxevals=maxevals)[1]
end