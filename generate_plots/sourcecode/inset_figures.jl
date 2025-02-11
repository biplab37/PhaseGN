using PhaseGN, PGFPlotsX, LaTeXStrings
using Plots

T = 1.0

param = Parameters(Λ=5.0, κ=kappa01(5.0))

q_list = [0.0, 0.5, 1.0, 1.2, 1.5, 2.0]

threshold_values = [sqrt(4 * mass_k(T, 0.0, q, param)^2 + q^2) for q in q_list]

values_near_threshold = [q-0.1:0.01:q+0.1 for q in q_list[2:end]]

ωrange = sort(union(10 .^ (-2:0.02:1.2), threshold_values, values_near_threshold...))

phases_q = zeros(length(q_list), length(ωrange))

function delta_phi(ω, q, T, param)
    m = mass_k(T, 0.0, q, param)
    impi = PhaseGN.imagpart_phi_q_refactored_m(ω, T, 0.0, q, m, param)
    repi = Π0_phi(T, 0.0, param) - PhaseGN.realpart(PhaseGN.imagpart_phi_q_refactored_m, ω, T, 0.0, q, m, param)
    return atan(impi, repi)
end

delta_phi(10.0, 0.0, 1.0, param)

function phase_shift_data(ωrange, q, T, param)
    phase_shift = delta_phi.(ωrange, q, T, param)
    return phase_shift
end

Threads.@threads for i in eachindex(q_list)
    phases_q[i, :] = phase_shift_data(ωrange, q_list[i], T, param)
end

# # save data
save_dir = "/home/biplab/Projects/PhaseGN/generate_plots/"
using DelimitedFiles
writedlm("$(save_dir)/data/phase_shift_data.dat", [["# T=0.85, ω, q=" q_list...]; ωrange phases_q'], ',')
using DelimitedFiles
data = readdlm("$(save_dir)/data/phase_shift_data.dat", ',', skipstart=1)

ωrange = data[:, 1]

delta_phi(4.0, 1.0, 0.1, param)

thr
boosted_phase_shift = [phase_shift_data(sqrt.(ωrange .^ 2 .+ q^2), 0.0, T, param) for q in q_list]

data_phase_shift = zeros(length(ωrange), length(q_list))

Threads.@threads for i in eachindex(ωrange)
    for j in eachindex(q_list)
        data_phase_shift[i, j] = last(delta_phi(ωrange[i], q_list[j], T, param))
    end
end
using Colors

colors = colormap("RdBu", length(q_list))

p = @pgf Axis(
    {
        xlabel = L"\omega^2",
        ylabel = L"\phi_{\varphi}(\omega, q)",
        xmode = "log",
        xmin = 0.05,
        xmax = 1.5e2,
        ymin = -0.1,
        ytick = [0, pi / 2, pi],
        yticklabels = ["0", L"\frac{\pi}{2}", L"\pi"],
    },
    Legend(
        ["q=$(q_list[i])M" for i in 1:4],
    ),
)

@pgf for i in [1, 2, 3, 4]
    lines = PlotInc(
        {
            no_marks,
            color = colors[i],
            style = "solid",
            # legend_entry = LaTeXString("q = $(q_list[i])")
        },
        Table(ωrange .^ 2, data[:, i+1])
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
    dashes = VLine(
        {
            style = "dashed",
            color = "black",
            opacity = 0.3
        },
        q_list[i]^2 + 4 * mass_k(T, 0.0, q_list[i], param)^2
    )
    push!(p, lines, dashes)
end

p

pgfsave("$(save_dir)/plots/5M/phase_shift_q_dependence.svg", p)

ωrange2 = sqrt(1.4):0.002:sqrt(2.4)

phases_q2 = zeros(length(q_list), length(ωrange2))

Threads.@threads for i in eachindex(q_list)
    phases_q2[i, :] = phase_shift_data(ωrange2, q_list[i], T, param)
end

p2 = @pgf Axis(
    {
    xlabel = L"\omega^2",
    width = "5cm",
    height = "5cm",
    # ylabel = L"\phi_{\varphi}(\omega, q)",
    # xmode = "log",
    xmin = 1.43,
    xmax = 2.3,
    ymin = -0.1,
    ytick = [0, pi / 2, pi],
    yticklabels = ["0", L"\frac{\pi}{2}", L"\pi"],
},
# Legend(
#     ["q=$(q_list[i])M" for i in 1:4],
# ),
)

@pgf for i in [1, 2, 3, 4]
    lines = PlotInc(
        {
            no_marks,
            # color = "black",
            style = "solid",
            # legend_entry = LaTeXString("q = $(q_list[i])")
        },
        Table(ωrange2 .^ 2, phases_q2[i, :])
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
    dashes = VLine(
        {
            style = "dashed",
            color = "black",
            opacity = 0.3
        },
        q_list[i]^2 + 4 * mass_k(T, 0.0, q_list[i], param)^2
    )
    push!(p2, lines, dashes)
end

p2


pgfsave("$(save_dir)/plots/5M/phase_shift_q_dependence_inset1.svg", p2)

# plot(ωrange, phases_q', xaxis=:log, lab=["q=0.0" "q=100M"], legend=:topleft, xlabel=L"\omega^2/M^2", ylabel=L"\phi_{\varphi}(\omega, q)")

# savefig("$(save_dir)/plots/phase_shift_q_dependence_tail.png")

## change in width linear scale

## Inset 2 near the end
ωrange3 = sqrt(99.5):0.002:sqrt(103)

phases_q3 = zeros(length(q_list), length(ωrange3))

Threads.@threads for i in eachindex(q_list)
    phases_q3[i, :] = phase_shift_data(ωrange3, q_list[i], T, param)
end

p3 = @pgf Axis(
    {
    xlabel = L"\omega^2",
    width = "5cm",
    height = "5cm",
    # ylabel = L"\phi_{\varphi}(\omega, q)",
    # xmode = "log",
    xmin = 100,
    xmax = 102.5,
    ymin = -0.05,
    # ymax=3.2,
    # ytick = [0, pi/2, pi],
    # yticklabels = ["0", L"\frac{\pi}{2}", L"\pi"],
},
# Legend(
#     ["q=$(q_list[i])M" for i in 1:4],
# ),
)

@pgf for i in [1, 2, 3, 4]
    lines = PlotInc(
        {
            no_marks,
            # color = "black",
            style = "solid",
            # legend_entry = LaTeXString("q = $(q_list[i])")
        },
        Table(ωrange3 .^ 2, phases_q3[i, :])
    )
    # boosted = PlotInc(
    #     {
    #         no_marks,
    #         color = "red",
    #         style = "dashed",
    #     # legend_entry = LaTeXString("q = $(q_list[i])")
    # },
    # Table(ωrange.^2, boosted_phase_shift[i])
    # )
    dashes = VLine(
        {
            style = "dotted",
            color = "black",
            opacity = 0.8
        },
        q_list[i]^2 + 4 * mass_k(T, 0.0, q_list[i], param)^2 + 4 * param.Λ^2
    )
    push!(p3, lines, dashes)
end

p3

pgfsave("$(save_dir)/plots/5M/phase_shift_q_dependence_inset2.svg", p3)

q_list[1]^2 + 4 * mass_k(T, 0.0, q_list[1], param)^2 + 4 * param.Λ^2
