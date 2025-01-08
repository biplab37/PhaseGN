using PhaseGN

param = Parameters(Λ=5.0, κ=0.045)

function delta_phi(ω, q, T, param)
    m = mass_k(T, 0.0, q, param)
    impi = PhaseGN.imagpart_phi_q_refactored_m(ω, T, 0.0, q, m, param)
    repi = - PhaseGN.realpart(PhaseGN.imagpart_phi_q_refactored_m, ω, T, 0.0, q, m, param)
    fullrepi = repi + Π0_phi(T, 0.0, param)
    ph_r = -atan(impi, -repi)
    ph_tot = atan(impi, fullrepi)
    ph_sc = ph_tot - ph_r
    return [ph_sc, ph_r, ph_tot]
end

function delta_sigma(ω, q, T, param)
    m = mass_k(T, 0.0, q, param)
    impi = PhaseGN.imagpart_sigma_q_refactored_m(ω, T, 0.0, q, m, param)
    repi = - PhaseGN.realpart(PhaseGN.imagpart_sigma_q_refactored_m, ω, T, 0.0, q, m, param)
    fullrepi = repi + Π0_sigma(T, 0.0, param)
    ph_r = -atan(impi, -repi)
    ph_tot = atan(impi, fullrepi)
    ph_sc = ph_tot - ph_r
    return [ph_sc, ph_r, ph_tot]
end

delta_phi(4.0, 1.0, 0.1, param)
delta_sigma(4.0, 1.0, 0.1, param)

ωrange = 10 .^(range(-2, stop=1.2, length=200))

t_list = [0.1, 0.7, 1.0, 2.0]
phases_k = zeros(length(ωrange), length(t_list), 3)

Threads.@threads for i in eachindex(ωrange)
    for j in eachindex(t_list)
        phases_k[i, j, :] = delta_phi(ωrange[i], 0.0, t_list[j], param)
    end
end

phases_k_sigma = zeros(length(ωrange), length(t_list), 3)

Threads.@threads for i in eachindex(ωrange)
    for j in eachindex(t_list)
        phases_k_sigma[i, j, :] = delta_sigma(ωrange[i], 0.0, t_list[j], param)
    end
end

@pgf gp = GroupPlot(
    {
        group_style = {
            group_size = "2 by 3",
            "x descriptions at=edge bottom",
            "y descriptions at=edge left",
            horizontal_sep = "0.6cm",
            vertical_sep = "0.4cm",
        },
        height = "4.2cm",
        width = "8.4cm",
        xmode = "log",
        xmin = 1e-3,
        # xmax = 1.5e2,
        xlabel = L"s = \omega^2",
        # ytick = raw"{-3.14, -1.57, 0, 1.57. 2.14}",
        # yticklabels=raw"{$-\pi$, $-\frac{\pi}{2}$, 0, $\frac{\pi}{2}$, $\pi$}",
        ytick = [-pi, -pi / 2, 0, pi / 2, pi],
        yticklabels = [L"-\pi", L"-\frac{\pi}{2}", "0", L"\frac{\pi}{2}", L"\pi"],
    },
    {
        ylabel = L"\phi_{\sigma, \rm{r}}",
        y_label_style = {at = "(axis description cs:-0.1,.5)",
            anchor = "south",
        },
    },
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,1,1])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,2,1])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,3,1])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,4,1])
    ),
    {
        ylabel = L"\phi_{\varphi, \rm{r}}",
        y_label_style = {at = "(axis description cs:1.01,.5)",
            anchor = "north",
        }
    },
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 1, 1])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 2, 1])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 3, 1])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 4, 1])
    ),
    {
        legend_pos = "south west",
        ylabel = L"\phi_{\sigma, \rm{sc}}",
        y_label_style = {at = "(axis description cs:-0.1,.5)",
            anchor = "south",
        }
    },
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,1,2])
    ),
    LegendEntry(L"T = 0.0"),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,2,2])
    ),
    LegendEntry(L"T = 0.7M"),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,3,2])
    ),
    LegendEntry(L"T = 1.0M"),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,4,2])
    ),
    LegendEntry(L"T = 2.0M"),
    {
        ylabel = L"\phi_{\varphi, \rm{sc}}",
        y_label_style = {at = "(axis description cs:1.01,.5)",
            anchor = "north",
        }
    },
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 1, 2])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 2, 2])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 3, 2])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 4, 2])
    ),
    {
        ylabel = L"\phi_{\sigma, \rm{tot}}",
        y_label_style = {at = "(axis description cs:-0.1,.5)",
            anchor = "south",
        },
        ymax = 3.3,
    },
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,1,3])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,2,3])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,3,3])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k_sigma[:,4,3])
    ),
    {
        ylabel = L"\phi_{\varphi, \rm{tot}}",
        y_label_style = {at = "(axis description cs:1.01,.5)",
            anchor = "north",
        }
    },
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 1, 3])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 2, 3])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 3, 3])
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases_k[:, 4, 3])
    ),
)

pgfsave("$(save_dir)/plots/5M/phases_q0.pdf", gp)