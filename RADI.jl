module RADI

"""
    timeprep(t_start, t_end, t_interval, saveperXsteps)

Prepare vectors of model timesteps and savepoints. All time units are in days.
"""
function timeprep(t_start, t_end, t_interval, saveperXsteps)
    timesteps = t_start:t_interval:t_end
    ntps = length(timesteps)
    savepoints = collect(1:saveperXsteps:ntps)
    # Save final timepoint if it's not already in the list
    if !(ntps in savepoints)
        append!(savepoints, ntps)
    end
    nsps = length(savepoints)
    return timesteps, savepoints, ntps, nsps
end # function timeprep

"""
    model()

Run the RADI model.
"""
function model(t_start, t_end, t_interval, saveperXsteps)

    timesteps, savepoints, ntps, nsps = timeprep(t_start, t_end, t_interval, saveperXsteps)
    sp = 1

    # Set water column conditions
    oxy_w = 400.0 # micromol/kg

    # Create arrays for modelled variables
    depths = LinRange(0.0, 20.0, 21) # in cm
    ndepths = length(depths)
    oxy = fill(oxy_w, (ndepths,))

    # Create arrays for saved variables
    oxy_save = fill(NaN, (ndepths, nsps))

    # Run RADI run: main model loop
    for (i, t) in enumerate(timesteps)

        # Model changes in variables
        oxy += depths .+ 0.1

        # Save output if we are at a savepoint
        if i in savepoints
            oxy_save[:, sp] = oxy
            sp += 1
        end # if
    end # for t

    return depths, oxy_save

end # function model

say_RADI() = println("RADI done!")

end # module RADI
