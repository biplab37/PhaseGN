using PhaseGN, PGFPlotsX, LaTeXStrings

function delta_phi(ω, q, T, param)
    impi = PhaseGN.imagpart_phi_q_refactored_k(ω, T, 0.0, q, param)
    repi = PhaseGN.realpart_3_k(PhaseGN.imagpart_phi_q_refactored_k, ω, T, 0.0, q, param)
    return atan(impi, repi)
end

function delta_integrand(ω, q, T, param)
    return delta_phi(ω, q, T, param) * q / ((exp(ω / T) - 1.0) * 2π^2)
end

function pressure_fluctuation(T, param)
    integrand(x) = delta_integrand(x[1], x[2], T, param) / T^3
    return PhaseGN.hcubature(integrand, [0.0, 0.0], [10 * T, 10 * T], abstol=1e-3)
end


using ProgressMeter
function pressure_fluctuations(trange, param)
    pressures = zeros(length(trange))
    @showprogress Threads.@threads for i in eachindex(trange)
        pressures[i] = pressure_fluctuation(trange[i], param)[1]
        println("Pressure fluctuation at T = ", trange[i], " is ", pressures[i])
    end
    return pressures
end

trange = range(0.1, 2.0, length=30)
param = Parameters()

@time pressures = pressure_fluctuations(trange, param)

using DelimitedFiles
writedlm("pressure_fluc_with_k2.dat", [["# full fluctuation pressure: T" "Pressure/T^3"]; trange pressures])
