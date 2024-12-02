## Function definitions for specific channels, their phases and the total phase.
#TODO: 

function phase_phi_q(temp, μ, ω, q, param::Parameters)
    return phase(imagpart_phi_q, Π0_phi_q, ω, temp, μ, q, param)
end

function phase_sigma_q(temp, μ, ω, q, param::Parameters)
    return phase(imagpart_sigma_q, Π0_sigma_q, ω, temp, μ, q, param)
end

function phasetot_phi_q(temp, μ, ω, q, param::Parameters)
    return phase_tot(imagpart_phi_q, Π0_phi_q, ω, temp, μ, q, param)
end

function phasetot_sigma_q(temp, μ, ω, q, param::Parameters)
    return phase_tot(imagpart_sigma_q, Π0_sigma_q, ω, temp, μ, q, param)
end

export phase_phi_q, phase_sigma_q, phasetot_phi_q, phasetot_sigma_q