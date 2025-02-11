using PhaseGN, PGFPlotsX, LaTeXStrings

function delta_sigma(ω, q, T, param)
    m = mass_k(T, 0.0, q, param)
    impi = PhaseGN.imagpart_sigma_q_refactored_m(ω, T, 0.0, q, m, param)
    repi = Π0_sigma(T, 0.0, param) - PhaseGN.realpart(PhaseGN.imagpart_sigma_q_refactored_m, ω, T, 0.0, q, m, param)
    return atan(impi, repi)
end

function mphi(temp, μ, param, bracket)
    m = mass_k(temp, μ, 0.0, param)
    func(ME) = Π0_phi(temp, μ, param) - PhaseGN.realpart(PhaseGN.imagpart_phi_q_refactored_m, ME, temp, μ, 0.0, m, param)
    return PhaseGN.bisection(func, bracket[1], bracket[2])
end

function Mphi_0(Λ, κ)
    param2 = Parameters(Λ=Λ, κ=κ)
    mphi(0.01, 0.0, param2, [0.0, 2.0])
end

function kappa01(Λ)
    f(κ) = Mphi_0(Λ, κ) - 0.1
    return PhaseGN.bisection(f, 0.0, 0.1)
end

function delta_integrand(ω, q, T, param)
    return delta_sigma(ω, q, T, param) * q / ((exp(ω / T) - 1.0) * 2π^2)
end

function pressure_fluctuation(T, param)
    integrand(x) = delta_integrand(x[1], x[2], T, param) / T^3
    return PhaseGN.hcubature(integrand, [0.0, 0.0], [min(10 * T, sqrt(5) * param.Λ), min(10 * T, param.Λ)], reltol=1e-1)
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