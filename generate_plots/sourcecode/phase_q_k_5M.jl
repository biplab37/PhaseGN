using PhaseGN, PGFPlotsX, LaTeXStrings

ωrange = 10 .^ (-2:0.01:1.2)

param = Parameters(Λ=5.0, κ=0.046)
μ = 0.0
q = 1.0

#Generate Data
tlist = 0.1:0.1:1.5
phases_phi = zeros(length(ωrange), length(tlist))

Threads.@threads for i in eachindex(tlist)
    phases_phi[:, i] = last.(delta_phi.(ωrange, q, tlist[i], param))
end

# save the data
using DelimitedFiles

writedlm("$(save_dir)/data/phase_q.dat", [ωrange phases_phi])

data = readdlm("$(save_dir)/data/phase_q.dat")

ωlist = data[:, 1]
ωbox = [[ωlist[1]]; ωlist; [ωlist[end]]]

axis = @pgf Axis(
    {
    set_layers,
    xmode = "log",
    view = "{-40}{25}",
    zmin = 0,
    zmax = 3.2,
    ymin = 0.2,
    ymax = 1.6,
    xmin = 8e-2,
    y_dir = "reverse",
    ytick = collect(tlist[3:2:end]),
    zticklabels = ["0", L"\pi/2", L"\pi"],
    ztick = [0, pi / 2, pi],
    xlabel = L"\omega/M",
    ylabel = L"T/M",
    zlabel = L"\phi_\varphi"
},
)

@pgf for i in 3:2:15
    curve = Plot3(
        {
            no_marks,
            style = {thick},
            color = "black"
        },
        Table(
            x=collect(ωlist),
            y=tlist[i] * ones(length(ωlist)),
            z=data[:, i+1]
        )
    )
    fill = Plot3(
        {
            draw = "none",
            fill = "black",
            fill_opacity = 0.15
        },
        Table(
            x=ωbox,
            y=tlist[i] * ones(length(ωbox)),
            z=[[0]; data[:, i+1]; [0]]
        )
    )
    push!(axis, curve, fill)
end

axis

pgfsave("$(save_dir)/plots/5M/phase_q.pdf", axis)