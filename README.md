# PhaseGN

[![codecov](https://codecov.io/gh/biplab37/PhaseGN/branch/main/graph/badge.svg?token=LZ1YO7D35E)](https://codecov.io/gh/biplab37/PhaseGN)

A Package to calculate the phase and realted quantities for $2+1$ dimensional Gross-Neveu Model.

## Installation
In the Julia REPL
```julia-repl
] add https://github.com/biplab37/PhaseGN.git
```
or using `Pkg`,
```julia
import Pkg; Pkg.add(url="https://github.com/biplab37/PhaseGN.git")
```
## Usage

```julia
using PhaseGN

p = Parameters()

phase_ϕ(temp,μ,ω,p)
```

For more detailed examples see the documentation.
