
using PGFPlotsX, PhaseGN
using LaTeXStrings

using DelimitedFiles

mott_data = readdlm("mott_momentas.dat")
trange = mott_data[:, 1]
mott_momentas = mott_data[:, 2]

function T_mott(mu, param)
    func(t) = M_phi(t, mu, param) - 2 * σ1(t, mu, param)
    return PhaseGN.bisection(func, 0.01, 1.0)
end
t_mott = T_mott(0.0, param)

p1 = @pgf Axis(
    {
        legend_pos = "north west",
        xlabel = "T/M",
        ylabel = L"q_{\rm Mott}/M",
        color = "black",
    },
    Plot(
        {
            color = "black",
            mark = "o",
            mark_options = "{scale=0.5}",
        },
        Table(trange, mott_momentas)
    ),
    # LegendEntry(L"q_{\rm Mott}"),
    # PlotInc(
    #     {no_marks, black},
    #     Table(trange, model.(trange))
    # ),
    # LegendEntry(L"\sqrt{2(T^2 - T_{\rm Mott}^2)}")
    )
pgfsave("$(save_dir)/plots/mott_momentas_new.pdf", p1)