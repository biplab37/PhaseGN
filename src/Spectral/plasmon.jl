## This file contains code related to the bulk and surface plasmons

function plasmon_frequency(T, mu, q, param)
    m = mass_k(T, mu, q, param)
    f(ω) = Π0_phi_m(T, mu, m, param) - realpart(imagpart_phi_q_refactored_m, ω, T, mu, q, m, param)
    return bisection(f, 0., 1.5)
end