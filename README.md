# Radi.jl

The one-dimensional reactive-advective-diffusive-irrigative diagenetic sediment module in Julia.

## Instructions

### Run the model in Julia

Clone or otherwise save this repo locally.

Next, prepare a setup script (e.g. [setup/IC_W29.jl](https://github.com/RADI-model/Radi.jl/blob/master/setup/IC_W29.jl)) that contains the initial conditions for the problem to be investigated.

Then, in Julia:

```julia
# Import the RADI model
include("Radi.jl")

# Run RADI for the first time
results = Radi.go("setup/IC_W29.jl")  # or a different setup script

# Access the different variables:
results[:savetimes]  # time of savepoints in years
results[:depths]  # model depths in m
results[:dO2]  # dissolved oxygen (rows = depths, columns = savetimes)
# and so on.

# Do a new run starting at the previous endpoint (overwrites previous results):
results = Radi.again(results)
# this uses internal settings from whichever file was last used with Radi.go().

# At any point, save the results as "results/IC_W29.mat":
Radi.save(results)  # overwrite original file, or...
Radi.save(results, "_more")  # ... append "_more" to the file name
# the "IC_W29" part of the filename is defined by `modelrun` in the setup file
```

### Plot the results in GNU Octave/MATLAB

Julia does not generate a plot of the results, but instead saves the results to a .mat file in the [`results/`](results) directory.  Use [`plot/everything.m`](plot/everything.m) to import and plot the results saved by Julia in GNU Octave/MATLAB.
