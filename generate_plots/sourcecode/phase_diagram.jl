using PhaseGN, PGFPlotsX, LaTeXStrings

trange = 0.01:0.001:1.1

param = Parameters(κ=0.0, Λ=5.0)

mus = PhaseGN.critical_line.(trange, param)

p = @pgf Axis({
        ylabel = L"T/M",
        xlabel = L"\mu/M",
        xmin = 0.0,
        # xmax = 1.1,
        ymin = 0.0,
        ymax = 0.8
    },
    Plot({
            no_marks,
        },
        Table(mus, trange)
    ),
    [raw"\node ", " at ", Coordinate(0.3, 0.3), "{\$\\bar{\\Phi}_1\\neq 0\$};"],
    [raw"\node ", " at ", Coordinate(0.9, 0.65), "{\$\\bar{\\Phi}_1 = 0\$};"],
    # [raw"\node ", " at ", Coordinate(0.8, 0.1), "{\$T_c=M/2\\ln 2\$};"],
    # VLine(
    #     {
    #         style="dashed",
    #         opacity=0.3,
    #     },
    #     1/2log(2)
    # )
    Plot(
        {
            only_marks,
        },
        Table([1.0], [0.0])
    )
)

pgfsave("critical_line.pdf", p)
