# Radi.jl

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.4680228-informational)](https://doi.org/10.5281/zenodo.4680228)

The one-dimensional reactive-advective-diffusive-irrigative diagenetic sediment module in Julia.

## Instructions

### Run the model in Julia

Clone or otherwise save this repo locally.  You need at least the main [Radi.jl](Radi.jl) file and the directory [modules](modules) plus all its contents.

Next, prepare a setup script (e.g. [setup/IC_W29.jl](setup/IC_W29.jl)) that contains the initial conditions for the problem to be investigated.

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

Julia does not generate a plot of the results, but instead saves the results to a .mat file in the [results](results) directory.  Use [plot/everything.m](plot/everything.m) to import and plot the results saved by Julia in GNU Octave/MATLAB.

## Citation

To cite this Julia implementation of RADI, please use:

> Humphreys, M. P., and Sulpis, O. (2021). **Radi.jl: the reactive-advective-diffusive-irrigative diagenetic sediment module in Julia.**  *Zenodo.*  [doi:10.5281/zenodo.4680228](https://doi.org/10.5281/zenodo.4680228).

[Bibtex version available here.](radi_jl_2021.bib)

The DOI above is the 'concept DOI' for all versions of Radi.jl; it will always automatically resolve to the latest release.  If you use a specific version of Radi.jl, please explicitly state the version, and switch to the relevant version-specific DOI:

  * v0.2: [doi.org/10.5281/zenodo.4680229](doi.org/10.5281/zenodo.4680229).

The rest of the citation remains the same.

There is also a manuscript in preparation that will describe RADI in detail.  Please check back here later for details.
