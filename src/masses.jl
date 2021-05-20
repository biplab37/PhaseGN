function σ1(temp,μ, κ=0.00, M=1.0)
    β = 1/temp
	
    σ₁(σ) = (1 + 2*cosh(β*μ)*exp(-β*σ) + exp(-2*β*σ) - exp(β*(M - σ + (π*κ/σ))))
	
    result = fzero(σ₁,1e-4,4)
	
    ## For the case when there is no zero in the interval
	if (4.0-result<1e-3 || result == 1e-4)
        result = 0.0
    end
	
    return result
end

export σ1