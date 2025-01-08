
using LaTeXStrings, PGFPlotsX
κ = kappa01(5.0)
p = Parameters(Λ=5.0, κ=κ)
trange = 0.01:0.01:1.15
μ = 0.0
m = σ1.(trange, μ, p);
Msigma = [mSigma(t, μ, p, [0.0, 3.0]) for t in trange];
Mphi = [mPhi(t, μ, p, [0.0, 3.0]) for t in trange];

p1 = @pgf Axis(
    {
        xlabel = "T/M",
        ylabel = "Mass [M]",
        xmin = 0.0,
        xmax = 1.15,
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
        Table(trange, maximum.(Mphis2))
    ),
    LegendEntry(L"M_\varphi"),
)

<<<<<<< Updated upstream
pgfsave("$(save_dir)/plots/5M/exciton_masses.pdf", p1)
=======
pgfsave("$(save_dir)/plots/5M/exciton_masses.pdf", p1)
>>>>>>> Stashed changes
