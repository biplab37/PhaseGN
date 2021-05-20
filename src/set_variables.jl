const Λ = [100.0]

function print_cutoff()
    println("The cutoff is "*string(Λ[1]))
end

function set_cutoff(new_Λ)
    Λ[1] = new_Λ
    println("The cutoff is set to "*string(Λ[1]))
end

export print_cutoff, set_cutoff