using PhaseGN, Plots, LaTeXStrings, PGFPlotsX, DelimitedFiles

function delta_phi(ω, q, T, param)
    impi = imagpart_phi_q(ω, T, 0.0, q, param)
    repi = Π0_phi(T, 0.0, param) - PhaseGN.realpart(imagpart_phi_q, ω, T, 0.0, q, param)
    return atan(impi, repi)
end

function delta_integrand(ω, q, T, param)
    return delta_phi(ω, q, T, param)*q/((exp(ω/T) -1.0)*2π^2)
end

function pressure_fluctuation(T, param)
    integrand(x) = delta_integrand(x[1], x[2], T, param)/T^3
    return PhaseGN.hcubature(integrand, [0.0, 0.0],[6*T, 6*T], reltol=1e-2)
end
using ProgressMeter
function pressure_fluctuations(trange, param)
    pressures = zeros(length(trange))
    @showprogress Threads.@threads for i in eachindex(trange)
        pressures[i] = pressure_fluctuation(trange[i], param)[1]
    end
    return pressures
end

trange = 0.1:0.02:2.0
param = Parameters()

@time pressures = pressure_fluctuations(trange, param)

plot(trange, pressures, marker=:circle)

writedlm("$(save_dir)/data/pressure_fluctuations.dat", [["# full fluctuation pressure: T" "Pressure/T^3"];trange pressures])
writedlm("pressure_fluctuations.dat", [["# full fluctuation pressure: T" "Pressure/T^3"];trange pressures])
function phase_zero(ω, T, mu, param)
    impi = imagpart_phi_q(ω, T, mu, 0.0, param)
    repi = Π0_phi(T, mu, param) - PhaseGN.realpart(imagpart_phi_q, ω, T, mu, 0.0, param)
    return atan(impi, repi)
end
pressure_boost = zeros(length(trange))
Threads.@threads for i in eachindex(trange)
    pressure_boost[i] = pressure(phase_zero, trange[i], 0.0, param)
end

p = @pgf Axis(
    {
        xlabel = L"T/M",
        ylabel = L"P/T^3",
        xmin = 0.0,
        ymin=0.0,
    },
    PlotInc(
        {
            color="black",
            # mark="*",
            no_marks,
        },
        Table(x=trange, y=pressures)
    ),
    LegendEntry("Full"),
    PlotInc(
        {
            no_marks,
            color="black",
            style="dashed",
        },
        Table(x=trange, y=pressure_boost)
    ),
    LegendEntry("Boosted")
)

pgfsave("$(save_dir)/plots/pressure_fluctuations.pdf", p)
