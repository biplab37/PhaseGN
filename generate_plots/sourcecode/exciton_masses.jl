using Plots
using LaTeXStrings, PGFPlotsX
param = Parameters(Λ=5.0, κ=kappa01(5.0))
trange = 0.9:0.001:1.0
μ = 0.0
m = mass_k.(trange, 0.0, 0.0, param)
Msigma = [mSigma(t, μ, param, [0.01, 3.0]) for t in trange];
Mphi = [mPhi(t, μ, p, [0.01, 2.0]) for t in trange];

Mphis2 = []
for i in 1:length(trange)
    mphis = mPhi_s2(trange[i], μ, param)
    push!(Mphis2, mphis)
end

Mphis2[1]

plot(trange,  [minimum.(Mphis2) maximum.(Mphis2)])

plot(trange, [2m, Msigma, Mphi2], label=["2m" L"M_\sigma" L"M_\varphi"])
vline!([0.98])
p = @pgf Axis(
    {
        xlabel = "T/M",
        ylabel = "Mass [M]",
        xmin = 0.0,
        xmax = 1.17,
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


plot(trange, [2m minimum.(Mphis2) maximum.(Mphis2)], xlims=(0.9,1))

pgfsave("exciton_masses.svg", p)