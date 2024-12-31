using PhaseGN, PGFPlotsX, LaTeXStrings

function delta_phi(ω, q, T, param)
    m = mass_k(T, 0.0, q, param)
    impi = PhaseGN.imagpart_phi_q_refactored_m(ω, T, 0.0, q, m, param)
    repi = Π0_phi(T, 0.0,param) -  PhaseGN.realpart(PhaseGN.imagpart_phi_q_refactored_m, ω, T, 0.0, q, m, param)
    return atan(impi, repi)
end

function mPhi(temp, μ, param, bracket)
    m = mass_k(temp, μ, 0.0, param)
    func(ME) = Π0_phi(temp, μ, param) - PhaseGN.realpart(PhaseGN.imagpart_phi_q_refactored_m, ME, temp, μ, 0.0, m, param)
    return PhaseGN.bisection(func, bracket[1], bracket[2])
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
    return delta_phi(ω, q, T, param)*q/((exp(ω/T) - 1.0)*2π^2)
end

function pressure_fluctuation(T, param)
    integrand(x) = delta_integrand(x[1], x[2], T, param)/T^3
    return PhaseGN.hcubature(integrand, [0.0, 0.0],[min(10*T, sqrt(5)*param.Λ), min(10*T, param.Λ)], reltol=1e-1)
end
param = Parameters(Λ=5.0, κ=kappa01(5.0))

using ProgressMeter


function pressure_fluctuations(trange, param, file)
    pressures = zeros(length(trange))
    @showprogress Threads.@threads for i in eachindex(trange)
        pressures[i] = pressure_fluctuation(trange[i], param)[1]
        writedlm(file, [trange[i] pressures[i]])
        println("Pressure fluctuation at T = ", trange[i], " is ", pressures[i])
    end
    return pressures
end

trange = range(0.1, 2.0, length=60)
param = Parameters(Λ=5.0, κ=kappa01(5.0))

file = open("pressure_fluc_with_k2_5M.dat", "a")

@time pressures = pressure_fluctuations(trange, param, file)
close(file)

using DelimitedFiles
writedlm("pressure_fluc_with_k2_1_5M.dat", [["# full fluctuation pressure: T" "Pressure/T^3"];trange pressures])

using Plots
plot(trange, pressures, label = "Pressure fluctuation", xlabel = "T", ylabel = "Pressure/T^3", title = "Pressure fluctuation with k = 1.5M", legend = :topleft)

data = readdlm("pressure_fluc_with_k2_1_5M.dat", skipstart = 1)

plot(data[:, 1], data[:, 2], marker=:circle, label = "", xlabel = "T", ylabel = "Pressure/T^3", title = "Pressure fluctuation", legend = :topleft)

function phase_zero(ω, T, mu, param)
    return delta_phi(ω, 0.0, T, param)
end

trange2 = range(0.01, 2.0, length=100)
pressure_boost = zeros(length(trange2))
Threads.@threads for i in eachindex(trange2)
    pressure_boost[i] = pressure2(phase_zero, trange2[i], 0.0, param)
end

pressure_boost2 = zeros(length(trange2))
Threads.@threads for i in eachindex(trange2)
    pressure_boost2[i] = boosted_fixed_pressure(phi_zero, trange2[i], 0.0, param)
end

function pressure2(func::Function, temp, μ, param::Parameters)
    integrand(s) = -(log(exp(sqrt(s) / temp) - 1) - sqrt(s) / temp) / (2 * π^2 * temp^2) * func(sqrt(s), temp, μ, param)
    int1(s) = integrand(1 / (1 - s)) / (1 - s)^2
    return PhaseGN.integrate(integrand, 0, 4*param.Λ^2)
end

plot(trange2, [pressure_boost pressure_boost2], label = "Pressure boost", xlabel = "T", ylabel = "Pressure/T^3", title = "Pressure fluctuation with k = 1.5M", legend = :topleft)



p1 = @pgf Axis(
    {
        xlabel = L"T",
        ylabel = L"\mathcal{P}_{\rm{fl, \varphi}}/T^3",
        legend_pos = "north east",
        xmin = 0.0, 
        xmax = 2.0,
        ymin = 0.0,
    },
    PlotInc(
        {
            no_marks,
            color = "black",
            style = "solid",
        },
        Table(x = data[:,1], y = data[:,2])
    ),
    LegendEntry("Full"),
    PlotInc(
        {
            no_marks,
            color = "black",
            style = "dashed",
        },
        Table(x = trange2, y = pressure_boost2)
    ),
    LegendEntry("Boosted")
)

pgfsave("pressure_comparison_5M_k.pdf", p1)

@time pressure_fluctuation(0.02, param)