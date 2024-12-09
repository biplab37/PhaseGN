using PhaseGN, PGFPlotsX, LaTeXStrings

q = 2.0
T = 0.85

param = Parameters()

ωrange = sort(union(10 .^ (-2:0.02:2.5), sqrt(4 * σ1(T, 0.0, param)^2 .+ q^2)))

phase_shift_dat = phase_shift_data(ωrange, q, T, param)

function boosted_phase_shifts(ωrange, q, T, param)
    return phase_shift_data(map(x->(x<=q) ? 0.0 : sqrt(x^2 - q^2), ωrange), 0.0, T, param)
end

delta_phi(0.0, 0.0, T, param)

ωrange2 = sqrt.(ωrange .^2 .+ q^2)
boosted_phase_shift_dat = boosted_phase_shifts(ωrange2, q, T, param)

plot(ωrange, phase_shift_dat, label="Phase shift", xaxis=:log)
plot!([0.01; ωrange2], [0.0; boosted_phase_shift_dat], label="Boosted phase shift", xaxis=:log)

comp = @pgf Axis(
    {
        xlabel = L"\omega^2",
        ylabel = L"\phi_{\varphi}(\omega, q)",
        xmode = "log",
        ytick = [0, pi/2, pi],
        yticklabels = ["0", L"\frac{\pi}{2}", L"\pi"],
        xmin = 0.005,
        xmax = 250.0^2,
        ymin = -0.1,
        ymax=3.3,
    },
    PlotInc(
        {
            no_marks,
            color = "black",
            style = "solid",
        },
        Table(x = ωrange.^2, y = phase_shift_dat)
    ),
    LegendEntry("Full"),
    PlotInc(
        {
            no_marks,
            color = "black",
            style = "dashed",
        },
        Table(x = [0.01; ωrange2.^2], y = [0.0; boosted_phase_shift_dat])
    ),
    LegendEntry("Boosted")
)

pgfsave("$(save_dir)/plots/boosted_phase_shift.pdf", comp)
## diff for a particular value of q
plot_diff = @pgf Axis(
    {
        xlabel = L"\omega^2",
        ylabel = L"\delta\phi_{\varphi}(\omega, q)",
        xmode = "log",
        # # xmin = -10.,
        # xmax = 10.0,
        ytick = [0, pi/2, pi],
        yticklabels = ["0", L"\frac{\pi}{2}", L"\pi"],
    },
    Plot(
        {
            no_marks,
            color = "black",
            style = "solid"
        },
        Table(x = ωrange .^2, y = phases_q[4, :] .- boosted_phase_shifts_q[4])

    )

)

pgfsave("$(save_dir)/plots/boosted_phase_difference.pdf", plot_diff)

## Difference between the two phase shifts at different q
q_list = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

ωrange = sort(union(10 .^ (-2:0.02:2.5), sqrt.(4 * σ1(T, 0.0, param)^2 .+ q_list.^2)))

phases_q = zeros(length(q_list), length(ωrange))

Threads.@threads for i in eachindex(q_list)
    phases_q[i, :] = phase_shift_data(ωrange, q_list[i], T, param)
end

boosted_phase_shifts_q = [boosted_phase_shifts(ωrange, q, T, param) for q in q_list]

p = @pgf Axis(
    {
        xlabel = L"\omega^2",
        ylabel = L"\phi_{\varphi}(\omega, q)",
        xmode = "log",
        # # xmin = -10.,
        # xmax = 10.0,
        ytick = [0, pi/2, pi],
        yticklabels = ["0", L"\frac{\pi}{2}", L"\pi"],
    },
    # Legend(
    #     ["q=$(q_list[i])M" for i in eachindex(q_list)],
    # )
)

@pgf for i in eachindex(q_list)
    lines = PlotInc(
        {
            no_marks,
            color = "black",
            style = "solid",
        },
        Table(x = (ωrange.^2), y = phases_q[i, :] .- boosted_phase_shifts_q[i])
    )
    push!(p, lines)
end

p
