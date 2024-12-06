
using LaTeXStrings
p = Parameters()
trange = 0.01:0.01:1
μ = 0.0
m = σ1.(trange, μ, p);
Msigma = [M_sigma(t, μ, p) for t in trange];
Mphi = [M_phi(t, μ, p) for t in trange];

p = @pgf Axis(
    {
        xlabel = "T/M",
        ylabel = "Mass [M]",
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
    },
    PlotInc(
        {
            no_marks,
            black,
            dashed,
        },
        Table(trange, 2 .* m)
    ),
    LegendEntry("2m"),
    PlotInc(
        {
            no_marks,
            black,
        },
        Table(trange, Msigma)
    ),
    LegendEntry(L"M_\sigma"),
    PlotInc(
        {
            no_marks,
            black,
            dashdotted,
        },
        Table(trange, Mphi)
    ),
    LegendEntry(L"M_\varphi"),
)

pgfsave("$(save_dir)/plots/exciton_masses.pdf", p)