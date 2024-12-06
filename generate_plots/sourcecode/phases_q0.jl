
ωrange = 10 .^ (-2:0.02:3);
phases = [phase_sigma(ω, 0.01, 0.0, Parameters()) for ω in ωrange];
phases_sc = first.(phases);
phases_tot = last.(phases);
phases_r = phases_tot - phases_sc;
phasep = [phase_phi(ω, 0.01, 0.0, Parameters()) for ω in ωrange];
phasep_sc = first.(phasep);
phasep_tot = last.(phasep);
phasep_r = phasep_tot - phasep_sc;
phases_1 = [phase_sigma(ω, 0.7, 0.0, Parameters()) for ω in ωrange];
phases1_sc = first.(phases_1);
phases1_tot = last.(phases_1);
phases1_r = phases1_tot - phases1_sc;
phasep1 = [phase_phi(ω, 0.7, 0.0, Parameters()) for ω in ωrange];
phasep1_sc = first.(phasep1);
phasep1_tot = last.(phasep1);
phasep1_r = phasep1_tot - phasep1_sc;
phases2 = [phase_sigma(ω, 1.5, 0.0, Parameters()) for ω in ωrange];
phases2_sc = first.(phases2);
phases2_tot = last.(phases2);
phases2_r = phases2_tot - phases2_sc;
phasep2 = [phase_phi(ω, 1.5, 0.0, Parameters()) for ω in ωrange];
phasep2_sc = first.(phasep2);
phasep2_tot = last.(phasep2);
phasep2_r = phasep2_tot - phasep2_sc;
phasep3 = [phase_phi(ω, 3.0, 0.0, Parameters()) for ω in ωrange];
phasep3_sc = first.(phasep3);
phasep3_tot = last.(phasep3);
phasep3_r = phasep3_tot - phasep3_sc;
phases3 = [phase_sigma(ω, 3.0, 0.0, Parameters()) for ω in ωrange];
phases3_sc = first.(phases3);
phases3_tot = last.(phases3);
phases3_r = phases3_tot - phases3_sc;
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
        xmin = 0.6e-4,
        xmax = 1.1e6,
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
        Table(ωrange .^ 2, phases_r)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases1_r)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases2_r)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases3_r)
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
        Table(ωrange .^ 2, phasep_r)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phasep1_r)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phasep2_r)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phasep3_r)
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
        Table(ωrange .^ 2, phases_sc)
    ),
    LegendEntry(L"T = 0.0"),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases1_sc)
    ),
    LegendEntry(L"T = 0.7M"),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases2_sc)
    ),
    LegendEntry(L"T = 1.5M"),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases3_sc)
    ),
    LegendEntry(L"T = 3.0M"),
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
        Table(ωrange .^ 2, phasep_sc)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phasep1_sc)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phasep2_sc)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phasep3_sc)
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
        Table(ωrange .^ 2, phases_tot)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases1_tot)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases2_tot)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phases3_tot)
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
        Table(ωrange .^ 2, phasep_tot)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phasep1_tot)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phasep2_tot)
    ),
    PlotInc(
        {
            no_marks,
        },
        Table(ωrange .^ 2, phasep3_tot)
    ),
)
pgfsave("$(save_dir)/plots/phases_q0.pdf", gp)