## Calculating exciton mass at small cutoffs,

using PhaseGN, Plots, LaTeXStrings, DelimitedFiles
using Roots

function mPhi(temp, μ, param, bracket)
    m = mass_k(temp, μ, 0.0, param)
    func(ME) = PhaseGN.realpart_3(PhaseGN.imagpart_phi_q_refactored_m, ME, temp, μ, 0.0, m, param)
    zero_list = find_zeros(func, bracket)
    if length(zero_list) == 0
        return 0.0
    else
        return min(zero_list...)
    end
end

function mSigma(temp, μ, param, bracket)
    m = σ1(temp, μ, param)
    func(ME) = PhaseGN.realpart_3(PhaseGN.imagpart_sigma_q_refactored_m, ME, temp, μ, 0.0, m, param)
    zero_list = find_zeros(func, bracket)
    if length(zero_list) == 0
        return 0.0
    else
        return zero_list[1]
    end
end

mPhi(0.1, 0.0, param, [0.0, 3.0])


## 

function MPhi_0(Λ, κ)
    param2 = Parameters(Λ=Λ, κ=κ)
    mPhi(0.01, 0.0, param2, [0.0, 2.0])
end


function kappa01(Λ)
    f(κ) = MPhi_0(Λ, κ) - 0.1
    return PhaseGN.bisection(f, 0.0, 0.1)
end

@time kappa01(5.0)

param = Parameters(Λ=5.0, κ=kappa01(5.0))
MPhi_0(5, kappa01(5.0))


function mPhi_q(temp, μ, q, param, bracket, m)
    func(ME) = PhaseGN.realpart_3(PhaseGN.imagpart_phi_q_refactored_m, ME, temp, μ, q, m, param)
    return PhaseGN.bisection(func, bracket[1], bracket[2])
end

mPhi_q(0.9, 0.0, 1.0, param, [0.0, 3.0], mass_k(0.9, 0.0, 1.0, param))

function mpq(T, q, param;)
    m = mass_k(T, 0.0, q, param)
    return mPhi_q(T, 0.0, q, param, [0.0, q + 1], m)
end

param = Parameters(Λ=5.0, κ=kappa01(5.0))
qrange = 0.1:0.1:5.0

mpqs = zeros(length(qrange))

Threads.@threads for i in eachindex(qrange)
    mpqs[i] = mpq(0.9, qrange[i], param)
end

mpqs .- qrange

thresolds = [sqrt(4 * mass_k(0.9, 0.0, q, param)^2 + q^2) for q in qrange]

plot(qrange, [-mpqs + thresolds], xlabel="q", ylabel="width", title="Width of bound state", lab="")
savefig("width_of_bound_state.pdf")

## MOtt momenta
function mmott_momenta(T, param)
    m(q) = mass_k(T, 0.0, q, param)
    thre_val(q) = sqrt(4 * m(q)^2 + q^2)
    func(q) = PhaseGN.Π0_phi_m(T, 0.0, m(q), param) - PhaseGN.realpart(PhaseGN.imagpart_phi_q_refactored_m, thre_val(q), T, 0.0, q, m(q), param)
    if func(0.0) < 0.0
        return 0.0
    end
    return PhaseGN.bisection(func, 0.0, 3.5)
end


mmott_momenta(0.983, param)

trange = union(0.85:0.05:0.95, 0.98:0.002:1.0, 1.01:0.05:1.3)
trange = union(0.85:0.05:1.1, 1.102:0.002:1.15, 1.152:0.05:1.4)


Moott_momenta = zeros(length(trange))
Threads.@threads for i in eachindex(trange)
    Moott_momenta[i] = mmott_momenta(trange[i], param)
end

plot(trange, Moott_momenta, marker=:circle, xlabel="T", ylabel="q", title="Mott momenta", lab="")
plot!(trange, t -> (t < 0.98) ? 0.0 : 10 * sqrt((t - 0.95)^2 - (0.98 - 0.95)^2))

using LsqFit

@. model(x, p) = p[1] * sqrt((x - 0.93)^2 - (0.98 - 0.93)^2)

fit = curve_fit(model, trange[5:end], Moott_momenta[5:end], [10, 0.98])

fit.param

savefig("mott_momenta.pdf")

plot!(0.98:0.01:1.5, x -> 10 * (x - 0.98), lab="")

mPhi(0.95, 0.0, param, [0.0, 2.0])
2 * mass_k(0.95, 0.0, 0.0, param)

fitcurve = map(x -> (x < 0.98 ? 0.0 : 5 * (x^1.5 - 0.98^1.5)^(2 / 3)), trange)
fitcurve2 = map(x -> (x < 0.98 ? 0.0 : 10 * sqrt((x - 0.95)^2 - (0.98 - 0.95)^2)), trange)
plot!(trange, fitcurve2, lab="fit curve")

using PGFPlotsX
p2 = @pgf Axis(
    {
        xlabel = L"T/M",
        ylabel = L"q_{\rm{Mott}}/M",
        xmin = 0.82,
        # xmax = 1.51,
        ymin = -0.2,
        legend_pos = "north west",},
    PlotInc(
        {
            # no_marks,
            only_marks,
            mark = "o",
            color = "black",
            style = "solid",
            mark_options = {scale = 0.75}
        },
        Table(x=trange, y=Moott_momenta)
    ),
    LegendEntry("Mott momenta"),
    PlotInc(
        {
            no_marks,
            # mark = "o",
            color = "black",
            style = "dashed",
        },
        Table(x=trange, y=fitcurve2)
    ),
    LegendEntry("Fit curve"),
)
pgfsave("mott_momenta.pdf", p2)

function mott_temperature(mu, param)
    func(t) = mPhi(t, mu, param, [0.0, 2.0]) - 2 * mass_k(t, mu, 0.0, param)
    return PhaseGN.bisection(func, 0.01, 1.0)
end

t_mott = mott_temperature(0.0, param)

mmott_momenta(0.976, param)

mPhi_q(0.93, 0.0, 0.0, param, [0.0, 2.0], mass_k(0.93, 0.0, 0.0, param)) - 2 * mass_k(0.93, 0.0, 0.0, param)

function delta_phi(ω, q, T, param)
    m = mass_k(T, 0.0, q, param)
    impi = PhaseGN.imagpart_phi_q_refactored_m(ω, T, 0.0, q, m, param)
    repi = Π0_phi(T, 0.0, param) - PhaseGN.realpart(PhaseGN.imagpart_phi_q_refactored_m, ω, T, 0.0, q, m, param)
    return atan(impi, repi)
end

function phaseshiftplot(q, T, param)
    thhresold = sqrt(4 * mass_k(T, 0.0, q, param)^2 + q^2)
    ωrange = sort(union(thhresold, range(0.5, 1.0, length=100)))
    phase = zeros(length(ωrange))
    Threads.@threads for i in eachindex(ωrange)
        phase[i] = delta_phi(ωrange[i], q, T, param)
    end
    plot(ωrange, phase, xlabel="ω", ylabel="δ", title="Phase shift", lab="")
    hline!([π])
    vline!([thhresold])
end

phaseshiftplot(0.0, 0.94, param)

function realpart_phi_k(temp, M, param)
    m = mass_k(temp, 0.0, 0.0, param)
    return PhaseGN.realpart_3(PhaseGN.imagpart_phi_q_refactored_m, M, temp, 0.0, 0.0, m, param)
end

function Mphi_plot(temp, μ, param)
    m = mass_k(temp, μ, 0.0, param)
    func(ME) = PhaseGN.realpart_3(PhaseGN.imagpart_phi_q_refactored_m, ME, temp, μ, 0.0, m, param)
    mrange = 0.0:0.001:3.0
    plot(mrange, func, lab="", xlims=(0.0, 3.0), ylims=(-0.15, 0.15), xlabel=L"\omega", fontfamily="Computer Modern", ylabel=L"Re$\Pi_\varphi(0, \omega)$")
    # vline!([2*m])
end

Mphi_plot(1.02, 0.0, param)
mPhi_s(0.93, 0.0, param, [0.0, 2.0])
2 * mass_k(0.92, 0.0, 0.0, param)

mPhi(0.92, 0.0, param, [0.0, 2.0])

kappa01(5.0)

function mPhi_2(temp, μ, param, guess)
    m = mass_k(temp, μ, 0.0, param)
    func(ME) = PhaseGN.realpart_3(PhaseGN.imagpart_phi_q_refactored_m, ME, temp, μ, 0.0, m, param)
    return find_zero(func, guess)
end

function mPhi_s2(temp, μ, param)
    m = mass_k(temp, μ, 0.0, param)
    func(ME) = PhaseGN.realpart_3(PhaseGN.imagpart_phi_q_refactored_m, ME, temp, μ, 0.0, m, param)
    return find_zeros(func, 0.1, 3.0)
end


mPhi_s(1.02, 0.0, param)

2 * mass_k(0.9805, 0.0, 0.0, param)

## create animation for t with Mphi_plot

scalefontsizes(1.3)

anim = @animate for t in 0.85:0.01:1.1
    Mphi_plot(t, 0.0, param)
    title!("T = $t M")
end

gif(anim, "Mphi_plot.gif", fps=2)