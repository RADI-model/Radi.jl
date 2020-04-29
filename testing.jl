module tst

using Colors, Plots, Profile, ProfileView
Plots.default(show=true)
Plots.closeall()

include("gsw_rho.jl")
include("RADI.jl")
import .RADI

# Define model timesteps (all times in years)
stoptime = 50.0 # how long to run for
interval = 5/128000 # duration of each model timestep
saveperXsteps = 128000 # save results in intervals of this many timesteps

# Import site-specific settings
include("IC_W29.jl")

function radiwrap(dO2_i, dtCO2_i, pfoc_i)
    @time depths, dO2, dtCO2, pfoc = RADI.model(stoptime, interval,
        saveperXsteps, z_max, z_res, dbl, phiInf, phi0, beta, lambda_b,
        lambda_i, T, S, P, dO2_w, dtCO2_w, dtPO4_w, Fpom, Fpom_r, Fpom_s,
        Fpom_f,rho_pom, dO2_i, dtCO2_i, pfoc_i)
    return depths, dO2, dtCO2, pfoc
end # function radiwrap

function radiplot(dO2_i, dtCO2_i, pfoc_i)
    depths, dO2, dtCO2, pfoc = radiwrap(dO2_i, dtCO2_i, pfoc_i)
    ntps = size(dO2)[2]
    cmap = colormap("RdBu", ntps)
    cs = ntps
    p1 = plot(depths*100, dO2[:, 1]*1e3, legend=false, c=cmap[cs])
    p2 = plot(depths*100, dtCO2[:, 1]*1e3, legend=false, c=cmap[cs])
    p3 = plot(depths*100, pfoc[:, 1], legend=false, c=cmap[cs])
    for sp in 2:ntps
        plot!(p1, depths*100, dO2[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p2, depths*100, dtCO2[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p3, depths*100, pfoc[:, sp], legend=false, c=cmap[cs-sp+1])
    end # for sp
    plot(p1, p2, p3, layout=(3, 1))
    return depths, dO2, dtCO2, pfoc
end # function radiplot

depths, dO2, dtCO2, pfoc = radiplot(dO2_i, dtCO2_i, pfoc_i)

# showprofile = false
# if showprofile
#     Profile.clear()
#     @profile depths, dO2, pfoc = RADI.model(RADIargs..., dO2_i, pfoc_i)
#     ProfileView.view()
#     RADI.say_RADI()
# else
    # depths, dO2, dtCO2, pfoc = radiplot(dO2_i, pfoc_i)
# end

end # module tst
