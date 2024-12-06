
using LaTeXStrings
temprange = 0.01:0.01:1.5
param = Parameters()
pres_s(t) = pressure(phasesc_phi, t, 0.0, param)
pres_r(t) = pressure(phaser_phi, t, 0.0, param)
pressures_s = [pres_s(t) for t in temprange];
pressures_r = [pres_r(t) for t in temprange];
p = @pgf Axis(
    {
        xlabel = L"T/M",
        ylabel = L"\mathcal{P}/T^3",
        xmin = 0.0,
    },
    PlotInc(
        {
            mark = "none",
            color = "black",
            dashed = "true",
        },
        Table(temprange, pressures_s)
    ),
    LegendEntry(L"\mathcal{P}_{\rm{\varphi,sc}}"),
    PlotInc(
        {
            mark = "none",
            color = "black",
            dashdotted = "true",
        },
        Table(temprange, pressures_r)
    ),
    LegendEntry(L"\mathcal{P}_{\rm{\varphi,r}}"),
    PlotInc(
        {
            mark = "none",
            color = "black",
        },
        Table(temprange, pressures_r .+ pressures_s)
    ),
    LegendEntry(L"\mathcal{P}_{\rm{\varphi,tot}}"),
    PlotInc(
        {
            mark = "none",
            color = "black",
            dotted = "true",
        },
        Table(temprange, pressures_poll)
    ),
    LegendEntry(L"\mathcal{P}_{\rm{pole}}"),
)

pgfsave("$(save_dir)/plots/pole_pressures.pdf", p)