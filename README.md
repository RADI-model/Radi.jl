# RADI.jl

1D Reaction-Advection-Diffusion-Irrigation (RADI) diagenetic sediment module in Julia.

## Instructions

Note that the current running instructions are still very much a temporary bodge job!

```julia
# This sets things up and runs the model from initial constant conditions:
include("RADI.jl")

# Now save the output variables in the main scope:
dO2 = radi.dO2;  # dissolved oxygen
dtCO2 = radi.dtCO2;  # dissolved inorganic carbon
pfoc = radi.pfoc;  # fast-degrading POC
psoc = radi.psoc;  # slow-degrading POC
proc = radi.proc;  # refractory POC

# With the command below you can do a new run starting at the previous endpoint:
depths, dO2, dtCO2, pfoc, psoc, proc = radi.profiles(
  dO2[:, end], dtCO2[:, end], pfoc[:, end], psoc[:, end], proc[:, end])

# Just repeat the exact same command again as many times as you like:
depths, dO2, dtCO2, pfoc, psoc, proc = radi.profiles(
  dO2[:, end], dtCO2[:, end], pfoc[:, end], psoc[:, end], proc[:, end])
```
