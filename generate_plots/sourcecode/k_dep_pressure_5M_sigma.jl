using PhaseGN, PGFPlotsX, LaTeXStrings

function delta_sigma(ω, q, T, param)
    m = mass_k(T, 0.0, q, param)
    impi = PhaseGN.imagpart_sigma_q_refactored_m(ω, T, 0.0, q, m, param)
    repi = Π0_sigma(T, 0.0, param) - PhaseGN.realpart(PhaseGN.imagpart_sigma_q_refactored_m, ω, T, 0.0, q, m, param)
    return atan(impi, repi)
end


function mPhi(temp, μ, param, bracket)
    m = mass_k(temp, μ, 0.0, param)
    func(ME) = Π0_phi(temp, μ, param) - PhaseGN.realpart(PhaseGN.imagpart_sigma_q_refactored_m, ME, temp, μ, 0.0, m, param)
    return PhaseGN.bisection(f  unc, bracket[1], bracket[2])
end

function MPhi_0(Λ, κ)
    param2 = Parameters(Λ=Λ, κ=κ)
    mPhi(0.01, 0.0, param2, [0.0, 2.0])
end

function kappa01(Λ)
    f(κ) = MPhi_0(Λ, κ) - 0.1
    return PhaseGN.bisection(f, 0.0, 0.1)
end

function delta_integrand(ω, q, T, param)
    return delta_sigma(ω, q, T, param) * q / ((exp(ω / T) - 1.0) * 2π^2)
end

function delta_integrand_gen(ω, q, T, param)
    del = delta_sigma(ω, q, T, param)
    return (del - sin(2*del)/2) * q / ((exp(ω / T) - 1.0) * 2π^2)
end

function pressure_fluctuation(T, param)
    integrand(x) = delta_integrand(x[1], x[2], T, param) / T^3
    return PhaseGN.hcubature(integrand, [0.0, 0.0], [min(10 * T, sqrt(5) * param.Λ), min(10 * T, param.Λ)], reltol=1e-1)
end

function pressure_fluctuation_gen(T, param)
    integrand(x) = delta_integrand_gen(x[1], x[2], T, param) / T^3
    return PhaseGN.hcubature(integrand, [0.0, 0.0], [min(10 * T, sqrt(5) * param.Λ), min(10 * T, param.Λ)], reltol=1e-1)
end
param = Parameters(Λ=5.0, κ=kappa01(5.0))

using ProgressMeter

@time pressure_fluctuation_gen(0.5, param)

function pressure_fluctuations(trange, param, file)
    pressures = zeros(length(trange))
    @showprogress Threads.@threads for i in eachindex(trange)
        pressures[i] = pressure_fluctuation(trange[i], param)[1]
        writedlm(file, [trange[i] pressures[i]])
        println("Pressure fluctuation at T = ", trange[i], " is ", pressures[i])
    end
    return pressures
end

function pressure_fluctuations_gen(trange, param, file)
    pressures = zeros(length(trange))
    @showprogress Threads.@threads for i in eachindex(trange)
        pressures[i] = pressure_fluctuation_gen(trange[i], param)[1]
        writedlm(file, [trange[i] pressures[i]])
        println("Pressure fluctuation at T = ", trange[i], " is ", pressures[i])
    end
    return pressures
end

trange = range(0.1, 2.0, length=60)
param = Parameters(Λ=5.0, κ=kappa01(5.0))

file = open("pressure_fluc_with_k2_5M_25dec3_sigma.dat", "a")
file2 = open("pressure_fluc_with_k2_5M_25dec3_gen_sigma.dat", "a")
using DelimitedFiles
using Plots

@time pressures = pressure_fluctuations(trange, param, file)
@time pressures_gen = pressure_fluctuations_gen(trange, param, file2)
plot(trange, [pressures, pressures_gen])
close(file)
close(file2)


file1 = open("pressure_fluc_with_k2_5M_25dec31_sigma.dat", "a")
file21 = open("pressure_fluc_with_k2_5M_25dec3_gen1_sigma.dat", "a")

trange2 = range(0.03, 0.09, length=5)
@time pressures_smallT = pressure_fluctuations(trange2, param, file1)
@time pressures_smallT_gen = pressure_fluctuations_gen(trange2, param, file21)
close(file1)
close(file21)

writedlm("pressure_fluc_with_k2_1_5M.dat", [["# full fluctuation pressure: T" "Pressure/T^3"]; trange pressures])


function phase_zero(ω, T, mu, param)
    return delta_sigma(ω, 0.0, T, param)[3]
end
function pressure2(func::Function, temp, μ, param::Parameters)
    integrand(s) = -(log(exp(sqrt(s) / temp) - 1) - sqrt(s) / temp) / (2 * π^2 * temp^2) * func(sqrt(s), temp, μ, param)
    int1(s) = integrand(1 / (1 - s)) / (1 - s)^2
    return PhaseGN.integrate(integrand, 0, 4 * param.Λ^2)
end
function boosted_fixed_pressure(phi, temp, mu, param)
    func(s) = (log((1 - exp(sqrt(param.Λ^2 + s) / temp)) / (1 - exp(sqrt(s) / temp))) * temp - (sqrt(param.Λ^2 + s) - sqrt(s))) * phi(sqrt(s), temp, mu, param) / (8 * π^2 * temp^3)
    return PhaseGN.integrate(func, 0.0, 5 * param.Λ^2)
end

trange3 = range(0.05, 2.0, length=60)
pressure_boost = zeros(length(trange3))
Threads.@threads for i in eachindex(trange3)
    pressure_boost[i] = pressure2(phase_zero, trange3[i], 0.0, param)
end

pressure_boost2 = zeros(length(trange3))
Threads.@threads for i in eachindex(trange3)
    pressure_boost2[i] = boosted_fixed_pressure(phase_zero, trange3[i], 0.0, param)
end

writedlm("fig_10_a_pressure_comparison_sigma.csv", pressure_boost2)

plot(trange3, [pressure_boost2], label="Pressure boost", xlabel="T", ylabel="Pressure/T^3", title="Pressure fluctuation with k = 1.5M", legend=:topleft)

writedlm("boosted_pressure_fluc_with_k2_5M.dat", [["# boosted fluctuation pressure: T" "Pressure/T^3"]; trange3 pressure_boost2])

data = vcat(pressures_smallT, pressures)
x_data = vcat(trange2, trange)
data_gen = vcat(pressures_smallT_gen, pressures_gen)

p2 = @pgf Axis(
    {
        xlabel = L"T/M",
        ylabel = L"\mathcal{P}_{\rm{fl, \sigma}}/T^3",
        legend_pos = "north east",
        xmin = 0.0,
        xmax = 2.0,
        ymin = 0.0,
    },
    PlotInc(
        {
            no_marks,
            color = "black",
            style = "dashed",
        },
        Table(x=x_data, y=data)
    ),
    LegendEntry("BU"),
    PlotInc(
        {
            no_marks,
            color = "black",
            style = "solid",
        },
        Table(x=x_data, y=data_gen)
    ),
    LegendEntry("gen. BU"),
        PlotInc(
        {
            no_marks,
            color = "black",
            style = "dotted",
        },
        Table(x=trange3, y=pressure_boost2)
    ),
    LegendEntry("Boosted")
)

pgfsave("k_dep_pressure_comparison_5M_k_gen_sigma.pdf", p2)
