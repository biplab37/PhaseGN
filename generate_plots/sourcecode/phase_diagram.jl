using PhaseGN, PGFPlotsX, LaTeXStrings

trange = 0.01:0.001:1.1

param = Parameters(κ=0.0)

mus = PhaseGN.critical_line.(trange, param)


p = @pgf Axis({
        xlabel = L"T/M",
        ylabel = L"\mu/M",
        xmin = 0.0,
        xmax = 1.1,
        ymin = 0.0,
    },
    Plot({
            no_marks,
        },
        Table(trange, mus)
    ),
    [raw"\node ", " at ", Coordinate(0.3, 0.5), "{\$\\bar{\\Phi}_1\\neq 0\$};"],
    [raw"\node ", " at ", Coordinate(0.9, 0.7), "{\$\\bar{\\Phi}_1 = 0\$};"]
)

pgfsave("$(save_dir)/plots/critical_line.pdf", p)
