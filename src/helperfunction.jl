"""
    bisection(func::Function,start::Number,finish::Number,iteration::Integer=20)
Finds the zero of a function `func` inside a given interval (`start`,`finish`). Might fail
if there are multiple zeros or no zeros inside the interval.
"""
function bisection(func::Function, start::Number, finish::Number, iteration::Integer=20)
    mid = (start + finish) / 2.0

    for _ in 1:iteration
        if func(mid) * func(finish) > 0
            finish = mid
        else
            start = mid
        end
        mid = (start + finish) / 2.0
    end
    return mid
end

"""
	DiracDelta(input::Float64,δ::Float64 = 1e-3)
Dirac Delta Function that can be used in a numerical integration.
"""

function DiracDelta(input, δ=0.05)
    return exp(-input^2 / (4 * δ)) / (2 * sqrt(pi * δ))
end

"""
	PrincipalValue(x::Float64,ϵ::Float64 = 1e-3)
Useful function in Principle Value Integration.
"""
function PrincipalValue(x, ϵ=1e-3)
    if abs(x) < ϵ
        return 0.0
    else
        return 1 / x
    end
end

"""
	numberF(temp::Float64, μ::Float64, Energy::Float64)
Returns the number density at a given `Energy`, chemical potential `μ` and `temp` for Fermionic Species.
"""
function numberF(temp, μ, Energy)
    β = 1 / temp
    return 1 / (1 + exp(β * (Energy - μ)))
end

# using QuadGK
function integrate(func::Function, start::Number, finish::Number; maxevals=10000)
    return hquadrature(func, start, finish, reltol=1e-3, maxevals=maxevals)[1]
    # return quadgk(func,start,finish,rtol=1e-3,maxevals=maxevals)[1]
end

function integrate(func, start::Vector, finish::Vector; maxevals=100000)
    return hcubature(func, start, finish, reltol=1e-3, maxevals=maxevals)[1]
end

function fzero(func::Function, guess)
    sol = nlsolve(x -> func(x...), [guess])
    # if sol.f_converged
    #     return sol.zero[1]
    # else
    #     @warn "No zero found. change the intial condition or increase number of iteration."
    # end
    return sol
end

function find_zero_interval_start(func::Function, a::Float64, b::Float64; tol::Float64=1e-6, max_iter::Int=100)
    """
    Find the start point 'c' of the interval (c,d) where func(x) = 0,
    given that func(x) >= 0 for all x in [a,b] and a < c <= d <= b.

    Parameters:
    func : Function
        The non-negative function to analyze.
    a, b : Float64
        The interval to search in.
    tol : Float64
        Tolerance for the result.
    max_iter : Int
        Maximum number of iterations.

    Returns:
    Float64
        The approximate value of 'c', or throws an error if not found.
    """

    is_zero(x, epsilon=1e-10) = abs(func(x)) < epsilon

    if is_zero(a)
        return a
    end

    if func(b) > π / 2
        b = 2b
        tol = 10 * tol
        @warn "Maximum value increased!"
    end

    # Check if the right endpoint is already zero
    if is_zero(b)
        # Search backwards from b to find where the function becomes non-zero
        left, right = a, b
        for _ in 1:max_iter
            if right - left < tol
                return right
            end
            mid = (left + right) / 2
            if is_zero(mid)
                right = mid
            else
                left = mid
            end
        end
        error("Maximum iterations reached while searching backwards")
    end

    # Original algorithm for when right endpoint is non-zero
    left, right = a, b

    for _ in 1:max_iter
        if right - left < tol
            return left
        end

        mid = (left + right) / 2
        f_left, f_mid = func(left), func(mid)

        if is_zero(left)
            return left
        elseif is_zero(mid)
            right = mid
        elseif f_left > 0 && f_mid > 0
            if f_left != f_mid
                guess = left - f_left * (right - left) / (f_mid - f_left)
                guess = max(left, min(right, guess))
            else
                guess = mid
            end

            if is_zero(guess)
                return guess
            elseif func(guess) > 0
                left = guess
            else
                right = guess
            end
        else
            error("Function appears to be negative in the interval")
        end
    end

    error("Maximum iterations reached without finding the zero interval")
end
