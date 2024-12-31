using PGFPlotsX, LaTeXStrings, PhaseGN

param = Parameters(Λ=5.0, κ=kappa01(5.0))
param0 = Parameters(Λ=5.0, κ=0.0)

srange = -2:0.01:2


T1 = 0.01
mu_list = [0.0, 0.9, 1.0, 1.1]


potentials = zeros(length(mu_list), length(srange))
potentials0 = zeros(length(mu_list), length(srange))

for (i, mu) in enumerate(mu_list)
    potentials[i, :] = PhaseGN.Omega_kappa.(srange, T1, mu, param)
end

for (i, mu) in enumerate(mu_list)
    potentials0[i, :] = PhaseGN.Omega_kappa.(srange, T1, mu, param0)
end

p = @pgf Axis(
    {
        xlabel = L"\bar{\Phi}_1/M",
        ylabel = L"\Omega/M^3",
        xmin = -1.8,
        xmax = 1.8,
        ymin = -0.12,
        ymax = 0.1,
        ytick = [-0.1, -0.05, 0.0, 0.05, 1.0],
        yticklabels = [L"- 0.1", L"- 0.05", L"0.0", L"0.05", L"1.0"],
    },
    PlotInc(
        {
            no_marks,
            black,
        },
        Table(srange, potentials[1, :])
    ),
    LegendEntry(L"\mu=0.0"),
    PlotInc(
        {
            no_marks,
            black,
            dashed,
        },
        Table(srange, potentials[2, :])
    ),
    LegendEntry(L"\mu=0.9M"),
    PlotInc(
        {
            no_marks,
            black,
            dashdotted,
        },
        Table(srange, potentials[3, :])
    ),
    LegendEntry(L"\mu=M"),
    PlotInc(
        {
            no_marks,
            black,
            dotted,
        },
        Table(srange, potentials[4, :])
    ),
    LegendEntry(L"\mu=1.1M"),
)

pgfsave("phase_transition1.pdf", p)

