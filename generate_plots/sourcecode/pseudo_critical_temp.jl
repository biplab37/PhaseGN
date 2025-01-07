using Plots, PhaseGN, UsefulFunctions, LaTeXStrings
# set_theme(:default)
# scalefontsizes(0.8)

using Roots
function second_derivative(f, x; dx=2e-2)
    return (f(x + dx) + f(x - dx) - 2 * f(x)) / dx^2
end
function sg_filter(f, x; dx=5e-2)
    return (-1 * f(x - 2 * dx) + 16 * f(x - dx) - 30 * f(x) + 16 * f(x + dx) - 1 * f(x + 2 * dx)) / (12 * dx^2)
end
function first_derivative(f, x; dx=1e-2)
    return (-2 * f(x - 2 * dx) - f(x + dx) + f(x + dx) + 2 * f(x + 2 * dx)) / (10 * dx)
end

function pseudo_critical_temp(mu, q, param)
    der2_mass(t, param) = sg_filter(x -> mass_k(x, mu, q, param), t)
    try
        return find_zero(x -> der2_mass(x, param), (0.1, 1.5))
    catch
        return 0.0
    end
end

param = Parameters(Λ=5.0, κ=0.046)
param2 = Parameters(Λ=5.0, κ=0.1)
param3 = Parameters(Λ=5.0, κ=1.0)

mu_list = 1.12:0.0001:1.13
mu_list2 = 0.0:0.02:1.25
trange = 0.01:0.005:1.2

temperatures = pseudo_critical_temp.(mu_list2, 0.0, param2)
temperatures2 = pseudo_critical_temp.(mu_list, 0.0, param2)
plot(temperatures, mu_list2, marker=:circle)
plot!(temperatures2, mu_list, marker=:circle)
length(mu_list2)

#Omega near the jump
mu1 = 1.12
mu2 = 1.13

trange2 = 0.3:0.05:0.8

srange = -1.5:0.01:1.5
plot()
for t in trange2
    plot!(srange, s -> PhaseGN.Omega(s, t, mu1, param2), label="")
end
plot!(xlabel="m")


## Chiral susceptibility 
# second derivative of the potential with respect to temperature
tlist = 0.3:0.002:0.6

function susceptibility(temp, mu, param)
    m(t) = mass_k(t, mu, 0.0, param)
    return sg_filter(m, temp, dx=5e-2)
end

sus = [susceptibility(t, 1.132, param2) for t in tlist]

plot(tlist, sus)

function all_zeros(mu, param)
    m(temp) = sg_filter(t -> mass_k(t, mu, 0.0, param), temp, dx=5e-2)
    return reduce_unique(find_zeros(m, 0.11, 0.6))
end

function reduce_unique(list)
    n = length(list)
    unique_list = [list[1]]
    if n == 1
        return unique_list
    end
    for i = 2:n
        if sum(abs.(unique_list .- list[i]) .< 1e-2) == 0
            append!(unique_list, list[i])
        end
    end
    return unique_list
end


murange3 = 1.12:0.0002:1.15
all_zeros(1.135, param2)

data_pseudo_critical = zeros(length(murange3), 3)

Threads.@threads for i in eachindex(murange3)
    temps = all_zeros(murange3[i], param2)
    n = length(temps)
    if n == 2
        data_pseudo_critical[i, 1] = temps[1]
        data_pseudo_critical[i, 2] = temps[2]
    elseif n == 1
        data_pseudo_critical[i, 1] = temps[1]
    elseif n == 3
        data_pseudo_critical[i, :] = temps
    else
        print(temps)
    end
end


scatter(murange3, data_pseudo_critical, marker=:circle)

# reshape the data

function reshape_data(data)
    reshaped_data_x = Float64[]
    reshaped_data_y = Float64[]
    for i in 1:size(data)[1]
        for x in data[i, :]
            if x != 0.0
                append!(reshaped_data_x, x)
                append!(reshaped_data_y, murange3[i])
            end
        end
    end
    return reshaped_data_x, reshaped_data_y
end

tr, mur = reshape_data(data_pseudo_critical)

perm = sortperm(tr)
mur2 = mur[perm]
tr2 = sort(tr)


plot(tr2, mur2, marker=:circle, fontfamily="Computer Modern", boxstyle=:border, lab="", xlabel=L"T/M", ylabel=L"\mu/M", grid=false)

mu_list5 = sort(union(1.13:0.002:1.25, 0.0:0.05:1.0, 1.0:0.01:1.12))
temperatures3 = pseudo_critical_temp.(mu_list5, 0.0, param2)

plot(temperatures3[1:end-13], mu_list5[1:end-13])
scatter!(tr2, mur2, marker=:circle)

tr3 = vcat(tr2, temperatures3[1:end-14])
mur3 = vcat(mur2, mu_list5[1:end-14])

perm2 = sortperm(tr3)
mur4 = mur3[perm2]
tr4 = sort(tr3)

plot([0.0; tr4], [(1 + sqrt(1 + 4 * π * 0.1)) / 2; mur4])

trange = union(0.1:0.01:0.6, 0.6:0.001:0.725)

mus = [PhaseGN.critical_line(t, param2) for t in trange]

plot!([0.0; trange], [1.0; mus])

using PGFPlotsX

p1 = @pgf Axis(
    {
        xlabel = L"T/M",
        ylabel = L"\mu/M",
        ymin = 0.0,
        xmin = 0.0,
    },
    PlotInc(
        {
            no_marks,
            color = "black",
            style = "solid",
        },
        Table([0.0; trange], [1.0; mus]),
    ),
    LegendEntry(L"\kappa=0"),
    PlotInc(
        {
            no_marks,
            color = "black",
            style = "dashed",
        },
        Table([0.0; tr4], [(1 + sqrt(1 + 4 * π * 0.1)) / 2; mur4]),
    ),
    LegendEntry(L"\kappa=0.1M^2"),
    # HLine(
    #     {
    #         style = "dotted",
    #         opacity = 0.5,
    #     },
    #     (1 + sqrt(1 + 4 * π * 0.1)) / 2
    # )
)

pgfsave("phase_diagram.svg", p1)

inset = @pgf Axis(
    {
        width = "5cm",
        height = "3cm"
    },
    Plot(
        {
        },
        Table(tr2, mur2)
    )
)

pgfsave("phase_diagram_inset.svg", inset)