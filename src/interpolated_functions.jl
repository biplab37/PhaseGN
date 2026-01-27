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

    return (x, y) -> (x > 2 * sqrt(param.Λ^2 + y^2 / 4 + param.M^2) || y > q_factor * param.Λ) ? 0.0 : imagpart_interp(x, y)
end

function imreparts_phi_interpolated(T, mu, param; o_length=100, q_length=100, q_factor=2.0)
    ω_max = 2 * sqrt(param.Λ^2 * (1 + q_factor^2 / 4) + param.M^2)
    orange = 2.0 * range(0.0, ω_max, length=o_length)
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
    repart(x, y) = (x > 2.0 * ω_max || y > q_factor * param.Λ) ? pi0_interp(y) : repart_interp(x, y)
    return impart, repart
end

function pressure_fl_phi_interpolated(T, mu, param)
    imp, rep = imreparts_phi_interpolated(T, mu, param, o_length=500, q_factor=1.0)

    del(om, q) = atan(imp(om, q), rep(om, q))

    pres_integrand(om, q) = q * numberB(T, mu, om) * del(om, q) / (2 * π^2)

    return integrate(x -> pres_integrand(x...), [0.0, 0.0], [sqrt(5) * param.Λ, param.Λ])
end




"""
    adaptive_grid(f, min_x, max_x; rtol=1e-4, atol=1e-4, max_depth=12)

Generates an adaptive grid for the function `f` over the interval `[min_x, max_x]`.
It recursively subdivides intervals where linear interpolation would likely fail to meet the error tolerances.
Returns sorted arrays `x_points` and `y_points`.
"""
function adaptive_grid(f, min_x, max_x; rtol=1e-4, atol=1e-4, max_depth=12)
    function get_points(x1, y1, x2, y2, depth)
        mid_x = (x1 + x2) / 2
        mid_y = f(mid_x)
        lin_y = (y1 + y2) / 2
        delta = abs(lin_y - mid_y)

        if depth < max_depth && (delta > atol && delta / max(abs(mid_y), 1.0) > rtol)
            left_seg = get_points(x1, y1, mid_x, mid_y, depth + 1)
            right_seg = get_points(mid_x, mid_y, x2, y2, depth + 1)
            return [left_seg..., right_seg...]
        else
            return [(x2, y2)]
        end
    end

    y_min = f(min_x)
    y_max = f(max_x)
    segments = get_points(min_x, y_min, max_x, y_max, 0)

    final_points = [(min_x, y_min), segments...]

    return map(first, final_points), map(last, final_points)
end

"""
    phase_shift_phi_interpolated_adaptive(T, mu, q, param)

Returns a function which returns phase shifts as a function of ω, using adaptive grid for interpolation.
"""
function phase_shift_phi_interpolated_adaptive_m(T, mu, q, m, param)
    ω_max = 2 * sqrt(param.Λ^2 + q^2 / 4 + m^2)

    # Adaptive grid for imaginary part at q
    orange, impart_data = adaptive_grid(om -> imagpart_phi_q_refactored_m(om, T, mu, q, m, param), 0.0, ω_max)
    impart_interp = linear_interpolation(orange, impart_data)
    impart_func(x) = (x > ω_max) ? 0.0 : impart_interp(x)

    # For realpart, we need imagpart at q=0 as well
    if q != 0.0
        m0 = mass_k(T, mu, 0.0, param)

        # Adaptive grid for imaginary part at q=0
        ω_max_0 = 2 * sqrt(param.Λ^2 + m0^2)
        orange0, impart_data0 = adaptive_grid(om -> imagpart_phi_q_refactored_m(om, T, mu, 0.0, m0, param), 0.0, ω_max_0)
        impart_interp0 = linear_interpolation(orange0, impart_data0)
        impart_func0(x) = (x > ω_max_0) ? 0.0 : impart_interp0(x)
    else
        impart_func0 = impart_func
    end

    # Fast wrapper for realpart integration
    function fast_imagpart(ν, _T, _mu, _q, _m, _param)
        if _q == q
            return impart_func(ν)
        elseif _q == 0.0
            return impart_func0(ν)
        else
            return imagpart_phi_q_refactored_m(ν, _T, _mu, _q, _m, _param)
        end
    end

    pi0 = Π0_phi_m(T, mu, m, param)

    orange2, repart_data = adaptive_grid(om -> pi0 - realpart(fast_imagpart, om, T, mu, q, m, param), 0.0, 2 * ω_max)
    repart_interp = linear_interpolation(orange2, repart_data)
    repart_func(x) = (x > 2 * ω_max) ? pi0 : repart_interp(x)

    return x -> atan(impart_func(x), repart_func(x))
end

function phase_shift_phi_interpolated_adaptive(T, mu, q, param)
    m = mass_k(T, mu, q, param)
    return phase_shift_phi_interpolated_adaptive_m(T, mu, q, m, param)
end

"""
    phase_shift_sigma_interpolated_adaptive(T, mu, q, param)

Returns a function which returns phase shifts as a function of ω, using adaptive grid.
"""
function phase_shift_sigma_interpolated_adaptive_m(T, mu, q, m, param)
    ω_max = 2 * sqrt(param.Λ^2 + q^2 / 4 + m^2)

    # Adaptive grid for imaginary part at q
    orange, impart_data = adaptive_grid(om -> imagpart_sigma_q_refactored_m(om, T, mu, q, m, param), 0.0, ω_max)
    impart_interp = linear_interpolation(orange, impart_data)
    impart_func(x) = (x > ω_max) ? 0.0 : impart_interp(x)

    # For realpart, we need imagpart at q=0 as well
    if q != 0.0
        m0 = mass_k(T, mu, 0.0, param)

        # Adaptive grid for imaginary part at q=0
        ω_max_0 = 2 * sqrt(param.Λ^2 + m0^2)
        orange0, impart_data0 = adaptive_grid(om -> imagpart_sigma_q_refactored_m(om, T, mu, 0.0, m0, param), 0.0, ω_max_0)
        impart_interp0 = linear_interpolation(orange0, impart_data0)
        impart_func0(x) = (x > ω_max_0) ? 0.0 : impart_interp0(x)
    else
        impart_func0 = impart_func
    end

    # Fast wrapper for realpart integration
    function fast_imagpart(ν, _T, _mu, _q, _m, _param)
        if _q == q
            return impart_func(ν)
        elseif _q == 0.0
            return impart_func0(ν)
        else
            return imagpart_sigma_q_refactored_m(ν, _T, _mu, _q, _m, _param)
        end
    end

    pi0 = Π0_sigma_m(T, mu, m, param)

    orange2, repart_data = adaptive_grid(om -> pi0 - realpart(fast_imagpart, om, T, mu, q, m, param), 0.0, 2 * ω_max)
    repart_interp = linear_interpolation(orange2, repart_data)
    repart_func(x) = (x > 2 * ω_max) ? pi0 : repart_interp(x)

    return x -> atan(impart_func(x), repart_func(x))
end

function phase_shift_sigma_interpolated_adaptive(T, mu, q, param)
    m = mass_k(T, mu, q, param)
    return phase_shift_sigma_interpolated_adaptive_m(T, mu, q, m, param)
end

"""
    imagpart_phi_interpolated_adaptive(T, mu, param; o_length=100, q_length=100, q_factor=2.0)

Returns interpolated function of the imaginary part of the polarization function using adaptive grid in ω.
"""
function imagpart_phi_interpolated_adaptive(T, mu, param; q_length=100, q_factor=2.0)
    ω_max = 2 * sqrt(2 * param.Λ^2 + param.M^2)
    qrange = range(0.0, q_factor * param.Λ, length=q_length)

    base_omegas = Float64[]
    # Use representative q values to build common grid
    for q_val in [0.0, qrange[div(end, 2)], qrange[end]]
        m_q = mass_k(T, mu, q_val, param)
        om_pts, _ = adaptive_grid(om -> imagpart_phi_q_refactored_m(om, T, mu, q_val, m_q, param), 0.0, ω_max)
        append!(base_omegas, om_pts)
    end

    common_orange = sort(unique(base_omegas))
    if common_orange[1] != 0.0
        pushfirst!(common_orange, 0.0)
    end
    if common_orange[end] != ω_max
        push!(common_orange, ω_max)
    end
    common_orange = sort(unique(common_orange))

    data_imagpart = zeros(length(common_orange), length(qrange))

    Threads.@threads for i in eachindex(qrange)
        m_q = mass_k(T, mu, qrange[i], param)
        data_imagpart[:, i] = [imagpart_phi_q_refactored_m(om, T, mu, qrange[i], m_q, param) for om in common_orange]
    end

    imagpart_interp = linear_interpolation((common_orange, qrange), data_imagpart)

    return (x, y) -> (x > 2 * sqrt(param.Λ^2 + y^2 / 4 + param.M^2) || y > q_factor * param.Λ) ? 0.0 : imagpart_interp(x, y)
end

function imreparts_phi_interpolated_adaptive(T, mu, param; q_length=100, q_factor=2.0)
    impart = imagpart_phi_interpolated_adaptive(T, mu, param, q_length=q_length, q_factor=q_factor)

    ω_max = 2 * sqrt(param.Λ^2 * (1 + q_factor^2 / 4) + param.M^2)
    qrange = range(0.0, q_factor * param.Λ, length=q_length)

    base_omegas = Float64[]
    for q_val in [0.0, qrange[div(end, 2)], qrange[end]]
        m_q = mass_k(T, mu, q_val, param)
        pi0 = Π0_phi_m(T, mu, m_q, param)
        om_pts, _ = adaptive_grid(om -> pi0 - realpart_with2args(impart, om, q_val, m_q, param), 0.0, 2 * ω_max)
        append!(base_omegas, om_pts)
    end

    common_orange_re = sort(unique(base_omegas))
    if common_orange_re[1] != 0.0
        pushfirst!(common_orange_re, 0.0)
    end
    if common_orange_re[end] != 2 * ω_max
        push!(common_orange_re, 2 * ω_max)
    end
    common_orange_re = sort(unique(common_orange_re))

    data_realpart = zeros(length(common_orange_re), length(qrange))
    data_pi0 = zeros(length(qrange))

    Threads.@threads for i in eachindex(qrange)
        m_q = mass_k(T, mu, qrange[i], param)
        pi0 = Π0_phi_m(T, mu, m_q, param)
        data_pi0[i] = pi0
        data_realpart[:, i] = [pi0 - realpart_with2args(impart, om, qrange[i], m_q, param) for om in common_orange_re]
    end

    pi0_interp = linear_interpolation(qrange, data_pi0)
    repart_interp = linear_interpolation((common_orange_re, qrange), data_realpart)
    repart(x, y) = (x > 2.0 * ω_max || y > q_factor * param.Λ) ? pi0_interp(y) : repart_interp(x, y)

    return impart, repart
end


function pressure_fl_phi_interpolated_adaptive(T, mu, param; kwargs...)
    imp, rep = imreparts_phi_interpolated_adaptive(T, mu, param; q_factor=1.0, kwargs...)

    del(om, q) = atan(imp(om, q), rep(om, q))

    pres_integrand(om, q) = q * numberB(T, mu, om) * del(om, q) / (2 * π^2)

    return integrate(x -> pres_integrand(x...), [0.0, 0.0], [sqrt(5) * param.Λ, param.Λ])
end


