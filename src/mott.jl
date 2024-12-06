## This file contains code related to Mott efect

function T_Mott(mu, param)
    func(t) = M_phi(t, mu, param) - 2 * σ1(t, mu, param)
    return bisection(func, 0.0, 2 * param.M)
end

function mu_Mott(temp, param)
    func(mu) = M_phi(temp, mu, param) - 2 * σ1(temp, mu, param)
    return bisection(func, 0.0, 2 * param.M)
end

function phase_at_thresold(q, temp, mu, param::Parameters)
    m = σ1(temp, mu, param)
    omega0 = sqrt(4m^2 + q^2)
    # points = [omega0 - 0.01, omega0, omega0+0.01]
    return π - (PhaseGN.phasetot_phi_q(temp, mu, omega0, q, param))
end

function mott_momenta(temp, mu, param; max=1.0)
    func(q) = phase_at_thresold(q, temp, mu, param)
    return find_zero_interval_start(func, 0.0, max)
end

export T_Mott, mu_Mott, mott_momenta
