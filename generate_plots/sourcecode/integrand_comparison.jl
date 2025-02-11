using PhaseGN, Plots

function delta_phi(ω, q, T, param)
    m = mass_k(T, 0.0, q, param)
    impi = PhaseGN.imagpart_phi_q_refactored_m(ω, T, 0.0, q, m, param)
    repi = Π0_phi(T, 0.0, param) - PhaseGN.realpart(PhaseGN.imagpart_phi_q_refactored_m, ω, T, 0.0, q, m, param)
    return atan(impi, repi)
end

function phase_shifts(ωrange, q, T, param)
    return delta_phi.(ωrange, q, T, param)
end

function boosted_phase_shifts(ωrange, q, T, param)
    ωrange2 = map(x -> (x <= q) ? 0.0 : sqrt(x^2 - q^2), ωrange)
    return delta_phi.(ωrange2, 0.0, T, param)
end

param = Parameters(Λ=5.0, κ=kappa01(5.0))

T1 = 0.5

q = 1.0
ωrange = sort(union(10 .^ (-1.5:0.02:1.2), sqrt(4 * mass_k(T1, 0.0, q, param)^2 + q^2)))

phase_shift_dat = zeros(length(ωrange))

Threads.@threads for i in eachindex(ωrange)
    phase_shift_dat[i] = delta_phi(ωrange[i], q, T1, param)
end

boosted_phase_shift_dat = boosted_phase_shifts(ωrange, q, T1, param)

integrand_full = phase_shift_dat .* q ./ ((exp.(ωrange ./ T1) .- 1.0) .* 2π^2)
integrand_boosted = boosted_phase_shift_dat .* q ./ ((exp.(ωrange ./ T1) .- 1.0) .* 2π^2)

plot(ωrange, [phase_shift_dat boosted_phase_shift_dat], label="Phase shift", xaxis=:log)

plot(ωrange, [integrand_full integrand_boosted], label=["full" "boosted"], xaxis=:log)


function boosted_fixed_pressure(phi, temp, mu, param)
    func(s) = (log((1 - exp(sqrt(param.Λ^2 + s) / temp)) / (1 - exp(sqrt(s) / temp))) * temp - (sqrt(param.Λ^2 + s) - sqrt(s))) * phi(sqrt(s), temp, mu, param) / (8 * π^2 * temp^3)
    return PhaseGN.integrate(func, 0.0, 5 * param.Λ^2)
end

phi_zero(ω, temp, mu, param) = delta_phi(ω, 0.0, temp, param)


boosted_phase_shift_dat2 = boosted_fixed_pressure(phi_zero, T1, 0.0, param)

T1, T2 = 0.5, 1.0

q = 1.0

ωrange = sort(union(10 .^ (-1.5:0.02:1.2), sqrt(4 * mass_k(T1, 0.0, q, param)^2 + q^2)))

phase_shift_dat1 = zeros(length(ωrange))
phase_shift_dat2 = zeros(length(ωrange))

Threads.@threads for i in eachindex(ωrange)
    phase_shift_dat1[i] = delta_phi(ωrange[i], q, T1, param)
    phase_shift_dat2[i] = delta_phi(ωrange[i], q, T2, param)
end
boosted_phase_shift_dat1 = boosted_phase_shifts(ωrange, q, T1, param)
boosted_phase_shift_dat2 = boosted_phase_shifts(ωrange, q, T2, param)

integrand_full1 = phase_shift_dat1 .* q ./ ((exp.(ωrange ./ T1) .- 1.0) .* 2π^2)
integrand_full2 = phase_shift_dat2 .* q ./ ((exp.(ωrange ./ T2) .- 1.0) .* 2π^2)

integrand_boosted1 = boosted_phase_shift_dat1 .* q ./ ((exp.(ωrange ./ T1) .- 1.0) .* 2π^2)
integrand_boosted2 = boosted_phase_shift_dat2 .* q ./ ((exp.(ωrange ./ T2) .- 1.0) .* 2π^2)

using PGFPlotsX, LaTeXStrings

gp1 = @pgf GroupPlot(
    {
        group_style = {
            group_size = "2 by 2",
            # "x descriptions at=edge bottom",
            # "y descriptions at=edge left",
        },
        xlabel = L"s = \omega^2-q^2",
    },
    Plot(
        {
            no_marks,
        },
        Table(x=ωrange .^ 2, y=phase_shift_dat1)
    ),
    Plot(
        {
            no_marks,
        },
        Table(x=ωrange .^ 2, y=phase_shift_dat1)
    ),
    Plot(
        {
            no_marks,
        },
        Table(x=ωrange .^ 2, y=phase_shift_dat1)
    ),
    Plot(
        {
            no_marks,
        },
        Table(x=ωrange .^ 2, y=phase_shift_dat1)
    )
)