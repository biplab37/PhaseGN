## Function definitions for specific channels, their phases and the total phase.
#TODO: 

function phase_phi(temp, μ, ω, q, param::Parameters)
    return phase(imagpart_phi_q, Π0_phi, ω, temp, μ, q, param)
end

function phase_sigma(temp, μ, ω, q, param::Parameters)
    return phase(imagpart_sigma_q, Π0_sigma, ω, temp, μ, q, param)
end

function phasetot_phi(temp, μ, ω, q, param::Parameters)
    return phase_tot(imagpart_phi_q, Π0_phi, ω, temp, μ, q, param)
end

function phasetot_sigma(temp, μ, ω, q, param::Parameters)
    return phase_tot(imagpart_sigma_q, Π0_sigma, ω, temp, μ, q, param)
end
