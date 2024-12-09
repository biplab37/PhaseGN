using PhaseGN, PGFPlotsX, LaTeXStrings
using Plots

T = 0.85

param = Parameters()

function delta_phi(ω, q, T, param)
    impi = imagpart_phi_q(ω, T, 0.0, q, param)
    repi = Π0_phi(T, 0.0, param) - PhaseGN.realpart(imagpart_phi_q, ω, T, 0.0, q, param)
    return atan(impi, repi)
end

function phase_shift_data(ωrange, q, T, param)
    phase_shift = delta_phi.(ωrange, q, T, param)
    return phase_shift
end

qmott = data[31,2]

q_list = [0.0, 1.5, 10.0]

threshold_values = [sqrt(4 * σ1(T, 0.0, param)^2 + q^2) for q in q_list]

ωrange = sort(union(10 .^(-2:0.01:2.5), threshold_values))

phases_q = zeros(length(q_list), length(ωrange))

Threads.@threads for i in eachindex(q_list)
    phases_q[i, :] = phase_shift_data(ωrange, q_list[i], T, param)
end

# save data
using DelimitedFiles
writedlm("$(save_dir)/data/phase_shift_data.dat", [["# T=0.85, ω, q=" q_list...]; ωrange phases_q'], ',')

data = readdlm("$(save_dir)/data/phase_shift_data.dat", ',', skipstart=1)

boosted_phase_shift = [phase_shift_data(sqrt.(ωrange .^2 .+ q^2), 0.0, T, param) for q in q_list]

p = @pgf Axis(
    {
        xlabel = L"\omega^2",
        ylabel = L"\phi_{\varphi}(\omega, q)",
        xmode = "log",
        xmin=0.005,
        xmax=250.0^2,
        ymin=-0.1,
        ytick = [0, pi/2, pi],
        yticklabels = ["0", L"\frac{\pi}{2}", L"\pi"],
    },
    Legend(
        ["q=$(q_list[i])M" for i in eachindex(q_list)],
    )
)

@pgf for i in eachindex(q_list)
    lines = PlotInc(
        {
            no_marks,
            color = "black",
            style = "solid",
            # legend_entry = LaTeXString("q = $(q_list[i])")
        },
        Table(data[:,1].^2, data[:,i+1])
    )
    # boosted = PlotInc(
    #     {
    #         no_marks,
    #         color = "red",
    #         style = "dashed",
    #         # legend_entry = LaTeXString("q = $(q_list[i])")
    #     },
    #     Table(ωrange.^2, boosted_phase_shift[i])
    # )
    push!(p, lines)
end

p

pgfsave("$(save_dir)/plots/phase_shift_q_dependence.pdf", p)