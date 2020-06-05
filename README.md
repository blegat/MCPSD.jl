# MCPSD

This is a translation of the [`mc_psd.m`](https://www.math.aau.at/or/Software/mc_psd.m) MATLAB function written by Franz Rendl.
A preliminary version of this MATLAB function also appears in Section 7 of [HRVW96].

[HRVW96] Helmberg, Christoph, Franz Rendl, Robert J. Vanderbei, and Henry Wolkowicz.
"*An interior-point method for semidefinite programming.*"
SIAM Journal on Optimization 6, no. 2 (1996): 342-361.

The solver can be used with JuMP as follows:
```julia
using PCPSD, JuMP
model = direct_model(MCPSD.Optimizer())
x = @variable(model, [1:10] in MCPSD.Elliptope(5))
@objective(model, Min, x[1] + x[3] + x[6] + x[7] + x[8] + x[10])
optimize!(model)
value.(x)
```
