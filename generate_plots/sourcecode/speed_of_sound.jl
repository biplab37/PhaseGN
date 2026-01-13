p = Parameters(κ=0.01, Λ=5.0)

μ = 0.0

pres(t) = pressure_MF(t, μ, p)
enrgy(t) = energy_MF(t, μ, p)

trange = 0.1:0.02:2.0

pressures = [pres(t) for t in trange]
energies = [enrgy(t) for t in trange]

dp = diff(pressures)
de = diff(energies)
cs2 = (dp ./ de)

plot(trange[1:end-1], cs2, marker=:circle)

p1 = @pgf Axis(
    {
        xlabel = L"T/M",
        ylabel = L"c_S^2",
        xmin = 0.0,
        xmax = 2.0,
        ymin = 0.0,
        ymax = 1.0,
    },
    Plot(
        {
            no_marks,
        },
        Table(trange[1:end-1], cs2)
    ),
    HLine({dashed}, 1 / 2)
)

using DelimitedFiles

writedlm("fig_3_speed_of_sound.csv", hcat(trange[1:end-1], cs2))

pgfsave("$(save_dir)/plots/5M/speed_of_sound.pdf", p1)
