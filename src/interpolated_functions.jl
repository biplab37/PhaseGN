## This file contains code with interpolated function for speedups.

"""
	phase_shift_phi_interpolated(T, mu, q, param)

Returns a function which returns phase shifts as a function of ω.
"""
function phase_shift_phi_interpolated_m(T, mu, q, m, param)
    ω_max = 2 * sqrt(param.Λ^2 + q^2 / 4 + m^2)
    orange = range(0, ω_max, length=200)
    impart_data = [imagpart_phi_q_refactored_m(om, T, mu, q, m, param) for om in orange]
    impart_interp = linear_interpolation(orange, impart_data)
    impart_func(x) = (x > ω_max) ? 0.0 : impart_interp(x)

    orange2 = 2 * orange
    pi0 = Π0_phi_m(T, mu, m, param)
    repart_data = [pi0 - realpart(imagpart_phi_q_refactored_m, om, T, mu, q, m, param) for om in orange2]
    repart_interp = linear_interpolation(orange2, repart_data)
    repart_func(x) = (x > 2 * ω_max) ? pi0 : repart_interp(x)

    return x -> atan(impart_func(x), repart_func(x))
end

function phase_shift_phi_interpolated(T, mu, q, param)
    m = mass_k(T, mu, q, param)
    return phase_shift_phi_interpolated_m(T, mu, q, m, param)
end

"""
	phase_shift_sigma_interpolated(T, mu, q, param)

Returns a function which returns phase shifts as a function of ω.
"""
function phase_shift_sigma_interpolated_m(T, mu, q, m, param)
    ω_max = 2 * sqrt(param.Λ^2 + q^2 / 4 + m^2)
    orange = range(0, ω_max, length=200)
    impart_data = [imagpart_sigma_q_refactored_m(om, T, mu, q, m, param) for om in orange]
    impart_interp = linear_interpolation(orange, impart_data)
    impart_func(x) = (x > ω_max) ? 0.0 : impart_interp(x)

    orange2 = 2 * orange
    pi0 = Π0_sigma_m(T, mu, m, param)
    repart_data = [pi0 - realpart(imagpart_sigma_q_refactored_m, om, T, mu, q, m, param) for om in orange2]
    repart_interp = linear_interpolation(orange2, repart_data)
    repart_func(x) = (x > 2 * ω_max) ? pi0 : repart_interp(x)

    return x -> atan(impart_func(x), repart_func(x))
end

function phase_shift_sigma_interpolated(T, mu, q, param)
    m = mass_k(T, mu, q, param)
    return phase_shift_sigma_interpolated_m(T, mu, q, m, param)
end

"""
	phase_shifts_phi_interpolated(T, mu, param)

Returns a two dimensional function of ω and q for the phase shifts of pseudo-scalar channel.
"""
function phase_shifts_phi_interpolated(T, mu, param; o_length=100, q_length=100, q_factor=2.0)
    ω_max = 2 * sqrt(2 * param.Λ^2 + param.M^2)

    orange = range(0.0, ω_max, length=o_length)
    qrange = range(0.0, q_factor * param.Λ, length=q_length)

    data_ps = zeros(length(orange), length(qrange))

    Threads.@threads for i in eachindex(qrange)
        func = phase_shift_phi_interpolated(T, mu, qrange[i], param)
        data_ps[:, i] = map(func, orange)
    end

    return linear_interpolation((orange, qrange), data_ps)
end

"""
	phase_shifts_sigma_interpolated(T, mu, param)

Returns a two dimensional function of ω and q for the phase shifts of scalar channel.
"""
function phase_shifts_sigma_interpolated(T, mu, param; o_length=100, q_length=100, q_factor=2.0)
    m = mass_k(T, mu, q, param)
    ω_max = 2 * sqrt(2 * param.Λ^2 + m^2)

    orange = range(0.0, ω_max, length=o_length)
    qrange = range(0.0, q_factor * param.Λ, length=q_length)

    data_ps = zeros(length(orange), length(qrange))

    Threads.@threads for i in eachindex(qrange)
        func = phase_shift_sigma_interpolated_m(T, mu, qrange[i], m, param)
        data_ps[:, i] = map(func, orange)
    end

    return linear_interpolation((orange, qrange), data_ps)
end

"""

Returns interpolated function of the imaginary part of the polarization function as a function of ω and q
"""
function imagpart_phi_interpolated(T, mu, param; o_length=100, q_length=100, q_factor=2.0)
    ω_max = 2 * sqrt(2 * param.Λ^2 + param.M^2)
    orange = range(0.0, ω_max, length=o_length)
    qrange = range(0.0, q_factor * param.Λ, length=q_length)

    data_imagpart = zeros(length(orange), length(qrange))

    Threads.@threads for i in eachindex(qrange)
        m_q = mass_k(T, mu, qrange[i], param)
        data_imagpart[:, i] = [imagpart_phi_q_refactored_m(om, T, mu, qrange[i], m_q, param) for om in orange]
    end

    imagpart_interp = linear_interpolation((orange, qrange), data_imagpart)

    return (x, y) -> (x>2*sqrt(param.Λ^2 + y^2/4 + param.M^2) || y>q_factor*param.Λ) ? 0.0 : imagpart_interp(x, y)
end

function imreparts_phi_interpolated(T, mu, param; o_length=100, q_length=100, q_factor=2.0)
    ω_max = 2 * sqrt(param.Λ^2*(1 + q_factor^2/4) + param.M^2)
    orange = 2.0*range(0.0, ω_max, length=o_length)
    qrange = range(0.0, q_factor * param.Λ, length=q_length)

    impart = imagpart_phi_interpolated(T, mu, param, o_length=o_length, q_length=q_length, q_factor=q_factor)
    data_realpart = zeros(length(orange), length(qrange))
    data_pi0 = zeros(length(qrange))

    Threads.@threads for i in eachindex(qrange)
        m_q = mass_k(T, mu, qrange[i], param)
        pi0 = Π0_phi_m(T, mu, m_q, param)
        data_pi0[i] = pi0
        data_realpart[:, i] = [pi0 - realpart_with2args(impart, om, qrange[i], m_q, param) for om in orange]
    end

    pi0_interp = linear_interpolation(qrange, data_pi0)

    repart_interp = linear_interpolation((orange, qrange), data_realpart)
    repart(x, y) = (x>2.0*ω_max || y>q_factor*param.Λ) ? pi0_interp(y) : repart_interp(x, y)
    return impart, repart
end

function pressure_fl_phi_interpolated(T, mu, param)
	imp, rep = imreparts_phi_interpolated(T, mu, param, o_length=500, q_factor=1.0)

	del(om, q) = atan(imp(om, q), rep(om, q))

	pres_integrand(om, q) = q*numberB(T, mu, om)*del(om, q)/(2*π^2)

	return integrate(x->pres_integrand(x...), [0.0, 0.0], [sqrt(5)*param.Λ, param.Λ])
end

