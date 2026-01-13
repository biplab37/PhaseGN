using PhaseGN, PGFPlotsX, LaTeXStrings
using Plots

T = 1.0

param = Parameters(Λ=5.0, κ=0.046)


qmott = data[31, 2]

q_list = [0.0, 0.5, 1.0, 1.2, 1.5, 2.0]

threshold_values = [sqrt(4 * mass_k(T, 0.0, q, param)^2 + q^2) for q in q_list]

values_near_threshold = [q-0.1:0.01:q+0.1 for q in q_list[2:end]]

ωrange = sort(union(10 .^ (-2:0.02:1.2), threshold_values, values_near_threshold...))

# phases_q = zeros(length(q_list), length(ωrange))

# Threads.@threads for i in eachindex(q_list)
#     phases_q[i, :] = phase_shift_data(ωrange, q_list[i], T, param)
# end

# # save data
# using DelimitedFiles
# writedlm("$(save_dir)/data/phase_shift_data.dat", [["# T=0.85, ω, q=" q_list...]; ωrange phases_q'], ',')
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

p = @pgf Axis(
    {
        xlabel = L"\omega^2",
        ylabel = L"\phi_{\varphi}(\omega, q)",
        xmode = "log",
        xmin = 0.05,
        xmax = 2e2,
        ymin = -0.1,
        ytick = [0, pi / 2, pi],
        yticklabels = ["0", L"\frac{\pi}{2}", L"\pi"],
    },
    Legend(
        ["q=$(q_list[i])M" for i in 1:4],
    ),
)

@pgf for i in 1:4
    lines = PlotInc(
        {
            no_marks,
            # color = "black",
            style = "solid",
            # legend_entry = LaTeXString("q = $(q_list[i])")
        },
        Table(ωrange .^ 2, data_phase_shift[:, i])
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


writedlm("fig_8_a_resurrection_bound_state.csv", hcat(ωrange.^2, data_phase_shift), ',')

pgfsave("$(save_dir)/plots/5M/phase_shift_q_dependence.pdf", p)

# plot(ωrange, phases_q', xaxis=:log, lab=["q=0.0" "q=100M"], legend=:topleft, xlabel=L"\omega^2/M^2", ylabel=L"\phi_{\varphi}(\omega, q)")

# savefig("$(save_dir)/plots/phase_shift_q_dependence_tail.png")

## change in width linear scale

