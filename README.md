# RADI.jl

1D Reaction-Advection-Diffusion-Irrigation (RADI) diagenetic sediment module in Julia.

## Instructions

Note that the current running instructions are still very much a temporary bodge job!

```julia
# This sets things up and runs the model from initial constant conditions:
julia> include("testing.jl")

# Now save the output variables in the main scope:
julia> oxy = tst.oxy
julia> poc = tst.poc

# With the command below you can do a new run starting at the previous endpoint:
julia> depths, oxy, poc = tst.radiplot(oxy[:, end], poc[:, end])

# Just repeat the exact same command again as many times as you like:
julia> depths, oxy, poc = tst.radiplot(oxy[:, end], poc[:, end])
```
