const Λ = [100.0]

"""
	print_cutoff()

Prints the momentum cutoff used in for regularisation
"""
function print_cutoff()
    println("The cutoff is "*string(Λ[1]))
end

"""
	set_cutoff(new_Λ)

Sets the momentum cutoff used in regularisation to be `new_Λ`
"""
function set_cutoff(new_Λ)
    Λ[1] = new_Λ
    println("The cutoff is set to "*string(Λ[1]))
end

export print_cutoff, set_cutoff
