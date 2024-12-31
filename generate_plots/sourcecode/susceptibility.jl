trange = union(0.1:0.1:0.6,0.62:0.001:0.8,0.8:0.1:1.5)

function second_derivative(f,x;dx=1e-3)
    return (f(x+dx) + f(x-dx) - 2f(x))/dx^2
end

function susceptibility(temp, param)
    m(t) = mass_k(t, 0.0, 0.0, param)
    return second_derivative(m, temp, dx=1e-4)
end

function mdt(temp, param)
    m(t) = mass_k(t, 0.0, 0.0, param)
    return UsefulFunctions._derivative(m, temp, h=1e-4)
end

m_dt = mdt.(trange, param2)
mdt2 = mdt.(trange, Parameters(Λ=5.0, κ=0.001))
sus = susceptibility.(trange, param2)

plot(trange, [m_dt])
plot(trange, [sus])
vline!([1/(2log(2))])
param

mass_k(0.1+0.01im, 0.0, 0.0, param2)

ForwardDiff.derivative(x->mass_k(x, 0.0, 0.0, param), 0.01)

using PGFPlotsX, LaTeXStrings

p3 = @pgf Axis(
    {
        xlabel=L"T/M",
        ylabel=L"\frac{dm}{dT}",
        xmin=0.1,
        xmax=1.,
    },
    Plot(
        {
            no_marks,
            color="black"
        },
        Table(trange, m_dt)
    ),
    VLine({
        style="dashed",
        opacity=0.5,
    },
    1/(2*log(2)))
)

pgfsave("susceptibility.pdf",p3)