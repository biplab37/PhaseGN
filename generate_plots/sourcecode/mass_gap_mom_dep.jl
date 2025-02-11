using PhaseGN, PGFPlotsX, LaTeXStrings

param = Parameters(Λ=5.0, κ=kappa01(5.0))
param2 = Parameters(Λ=5.0, κ=0.0)

trange = 0.01:0.01:2.0

k_list = [0.0, 2.5, 5.0]

masses_k = zeros(length(k_list), length(trange))
masses_k0 = zeros(length(k_list), length(trange))

Threads.@threads for i in eachindex(k_list)
    for (j, t) in enumerate(trange)
        # guess = try 
        #     masses_k[i, j-1]
        # catch
        #     0.0
        # end
        masses_k[i, j] = mass_k(t, 0.0, k_list[i], param)
        masses_k0[i, j] = mass_k(t, 0.0, k_list[i], param2)
    end
end

p = @pgf Axis(
    {
        xlabel = L"T/M",
        ylabel = L"m/M",
        xmin = 0.0,
        xmax = 2.0,
        ymin = 0.0,
    },
    PlotInc(
        {
            color = "black",
            # mark="*",
            no_marks,
            style = "dotted"
        },
        Table(x=trange, y=masses_k0[1, :])
    ),
    LegendEntry(L"m_0=0"),
    PlotInc(
        {
            color = "black",
            # mark="*",
            no_marks,
        },
        Table(x=trange, y=masses_k[1, :])
    ),
    LegendEntry("q=0"),
    PlotInc(
        {
            color = "black",
            # mark="*",
            no_marks,
            style = "dashed",
        },
        Table(x=trange, y=masses_k[2, :])
    ),
    LegendEntry("q=2.5M"),
    PlotInc(
        {
            color = "black",
            # mark="*",
            no_marks,
            style = "dashdotted",
        },
        Table(x=trange, y=masses_k[3, :])
    ),
    LegendEntry("q=5M"),)

1 / (2log(2))

pgfsave("mass_gap_momentum_dep.svg", p)

