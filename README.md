# PhaseGN

A Package to calculate the phase and realted quantities for $2+1$ dimensional Gross-Neveu Model.

## Installation
```julia-repl
] add https://github.com/biplab37/PhaseGN.git
```
## Usage

```julia
using PhaseGN

Λ = 2. # GeV
set_cutoff(Λ)

phase_ϕ(temp,μ,ω,κ)
```