using PhaseGN, PGFPlotsX

save_dir = "/home/biplab/Documents/Projects/PhaseGN/paper/generate_plots/"

# include all the .jl file in the dir sourcecode 
for file in readdir(save_dir * "sourcecode")
    if occursin(".jl", file)
        include("$(save_dir * "sourcecode")/$file")
    end
end