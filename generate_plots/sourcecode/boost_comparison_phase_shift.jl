T1 = 0.5

orange = 10 .^(-1.5:0.01:1.25)

function phase_shifts(ωrange, q, T, param)
    return delta_phi.(ωrange, q, T, param)
end

function boosted_phase_shifts(ωrange, q, T, param)
    ωrange2 = map(x -> (x <= q) ? 0.0 : sqrt(x^2 - q^2), ωrange)
    return delta_phi.(ωrange2, 0.0, T, param)
end

phases1 = phase_shifts(orange, 1.0, T1, param)

phases1_boosted = boosted_phase_shifts(orange, 1.0, T1, param)

T2 = 1.0
phases2 = phase_shifts(orange, 1.0, T2, param)

phases2_boosted = boosted_phase_shifts(orange, 1.0, T2, param)

plot(orange, [phases1, phases1_boosted], xaxis=:log)
plot(orange, [phases2, phases2_boosted], xaxis=:log)

writedlm("fig_9_phase_shift_diff.csv", hcat(orange, phases1, phases1_boosted, phases2, phases2_boosted))

function integrand_mes(ω, q, T, param)
    return (1.0/ (exp(ω / T) - 1.0)) * delta_phi(ω, q, T, param) / (4π^2)
end

function integrand_mes_boosted(ω, q, T, param)
	om = (ω<=q) ? 0.0 : sqrt(ω^2 - q^2)
    return (1.0/ (exp(ω / T) - 1.0)) * delta_phi(om, 0.0, T, param) / (4π^2)
end

ints1 = integrand_mes.(orange, 1.0, 0.5, param)
ints1_boosted = integrand_mes_boosted.(orange, 1.0, 0.5, param)

ints2 = integrand_mes.(orange, 1.0, 1.0, param)
ints2_boosted = integrand_mes_boosted.(orange, 1.0, 1.0, param)
plot(orange, [ints1 ints1_boosted], xaxis=:log)
plot(orange, [ints2 ints2_boosted], xaxis=:log)

writedlm("fig_9_integrand_diff.csv", hcat(orange, ints1, ints1_boosted, ints2, ints2_boosted))