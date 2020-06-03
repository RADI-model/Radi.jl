# Radi.jl

The one-dimensional reactive-advective-diffusive-irrigative diagenetic sediment module in Julia.

## Instructions

Still a work in progress, but getting there.

### Run the model in Julia

```julia
# This sets things up and runs the model from initial constant conditions:
include("Radi.jl")

# Now save the results dictionary in the main scope, for convenience:
results = Radi.results;

# Access the different variables:
results[:savetimes]  # time of savepoints in years
results[:depths]  # model depths in m
results[:dO2]  # dissolved oxygen (rows = depths, columns = savetimes)
# and so on.

# Do a new run starting at the previous endpoint (overwrites previous dictionary):
results = Radi.again(results)

# Save the new set of results as a .mat file in /results/:
Radi.save(results)  # overwrite original file, or...
Radi.save(results, "_more")  # ... append "_more" to the file name
```

### Plot the results in Octave/MATLAB

Use [`plot/everything.m`](plot/everything.m) to import and plot the results saved by Julia.
