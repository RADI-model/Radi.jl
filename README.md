# Radi.jl

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.4680228-informational)](https://doi.org/10.5281/zenodo.4680228)

The one-dimensional reactive-advective-diffusive-irrigative diagenetic sediment module in Julia.

- [Radi.jl](#radijl)
  - [Instructions](#instructions)
    - [Run the  single model in Julia](#run-the-single-model-in-julia)
    - [Run the ensemble model in Julia](#run-the-ensemble-model-in-julia)
    - [Run the optimization model in Julia](#run-the-optimization-model-in-julia)
  - [Citation](#citation)

## Instructions

Clone or otherwise save this repo locally.  You need at least the version of the runfile that you would like to run (e.g. runfile, runfile_ensemble or runfile_opt) file and the directory [modules](modules) plus all its contents.

Next, prepare a setup script (e.g. [setup/IC_W29.jl](setup/IC_W29.jl)) that contains the initial conditions for the problem to be investigated.

Next, please make sure the following packages are installed in julia (can be installed in jupyter notebook with the command: import Pkg; Pkg.add("PackageName")
or in the julia terminal by: entering the package manager by pressing "]" and pkg> add "PackageName")

The following packages are needed to run RADIv2:
Distributed, MAT, DifferentialEquations, BenchmarkTools, ProgressLogging.
The packages: Plots and JLD2 are also handy to plot output (Plots) and save RADIv2 outputs locally (JLD2)

RADIv2 uses Julia's DifferentialEquations.jl package to solve the transport equations. https://docs.sciml.ai/DiffEqDocs/stable/ provides a comprehensive overview of this package, including how to use it and its options.

### Run the single model in Julia

To run the single model, ensure that the correct initial conditions (IC) are specified in the `runfile`. If the required packages (see installation instructions above) are properly installed, the model should run without modification.  
To customize the simulation, you can modify the following line:


```julia
@time single_sol = solve(prob, Rosenbrock23(autodiff=false), save_everystep=false)
```

This line solves the system using the Rosenbrock23 solver. You may replace it with other solvers, adjust tolerances, or set additional options. For more details on customizing the solver and available methods, refer to the [DifferentialEquations.jl documentation](https://diffeq.sciml.ai/).

### Run the ensemble model in Julia

To run the ensemble version of the model, ensure that the correct IC are set in the runfile_ensemble. In the provided example script, the model perturbs parameters such as T, U, omegaCa, FPOM, FPIC, kfast, kslow, and the fractions of fast- and slow-degrading POM (Fpom_f and Fpom_s).

You can adapt the model to vary other parameters by passing them into the calculate_constants() function. This allows each trajectory in the ensemble to use its own set of constants. These values can be either randomly generated (as in the example) or manually specified.

The number of ensemble members is set with:

```julia
trajectories = x
```

You can tune the batch size for parallel execution depending on your system’s capacity. The batch mechanism helps avoid overloading system resources when using multithreading.

To enable multithreading, you must manually set the number of threads for Julia, as it defaults to one. This can be done as follows:

export JULIA_NUM_THREADS=8 in macOS and Linux, and JULIA_NUM_THREADS=8 in Windows. 

Finally, ensure that the solver is called with threading enabled:

```julia
batch_sols = solve(
            local_prob,
            Rosenbrock23(autodiff=false),
            EnsembleThreads(), #enables multithreading
            trajectories=this_batch_size,
            callback=callback,  
            save_everystep=false
        )
```

For more details on parallel ensemble simulations and solver customization, consult the [DifferentialEquations.jl documentation](https://diffeq.sciml.ai/) documentation.


### Run the optimization model

To run the optimization model, make sure the following Julia packages are installed: BlackBoxOptim and Optim. The optimization process relies on a custom-defined loss function that assigns penalties when the model fails to meet specified criteria. This setup uses an ensemble solver, allowing the loss to be evaluated over multiple trajectories. In the example script provided, penalties are applied if model outputs such as DIC, TA, and O₂ flux fall outside predefined target ranges. However, the loss function is fully customizable—users can define it to optimize any input or output variable of the model, depending on their objectives.

## Citation

To cite this Julia implementation of RADI, please use:

> Humphreys, M. P., and Sulpis, O. (2021). **Radi.jl: the reactive-advective-diffusive-irrigative diagenetic sediment module in Julia.**  *Zenodo.*  [doi:10.5281/zenodo.4680228](https://doi.org/10.5281/zenodo.4680228).

[Bibtex version available here.](radi_jl_2021.bib)

The DOI above is the 'concept DOI' for all versions of Radi.jl; it will always automatically resolve to the latest release.  If you use a specific version of Radi.jl, please explicitly state the version, and switch to the relevant version-specific DOI:

  * v0.3: [doi:10.5281/zenodo.5005650](https://doi.org/10.5281/zenodo.5005650).
  * v0.2: [doi:10.5281/zenodo.4680229](https://doi.org/10.5281/zenodo.4680229).

The rest of the citation remains the same.

There is also a manuscript in preparation that will describe RADI in detail.  Radi.jl v1.0 will be released following peer review.  Please check back here later for details. 