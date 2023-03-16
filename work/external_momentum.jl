using PhaseGN

p = Parameters(κ = 0.01)

# calculate phase with non zero external momentum
phaser(ω,q;T=0.01,μ=0.0,param=p) = phaser_ϕ_q(T,μ,ω,q,param)
phasetot(ω,q;T=0.01,μ=0.0,param=p) = phase_ϕ_q(T,μ,ω,q,param)[3]

using Plots
plot!(10 .^(-2:0.02:2.2),x->phasetot(x,0.2,T=0.7),xaxis=:log)