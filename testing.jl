module tst

using Colors, Plots
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

"Convenient wrapper function for running RADI."
function radiwrap(dO2_i, dtCO2_i, pfoc_i, psoc_i, proc_i)
    @time depths, dO2, dtCO2, pfoc, psoc, pfoc = RADI.model(stoptime, interval,
        saveperXsteps, z_max, z_res, dbl, phiInf, phi0, beta, lambda_b,
        lambda_f, lambda_s, lambda_i, T, S, P, dO2_w, dtCO2_w, dtPO4_w, Fpom,
        Fpom_r, Fpom_s, Fpom_f, rho_pom, dO2_i, dtCO2_i, pfoc_i, psoc_i, proc_i)
end # function radiwrap

"Run RADI and plot the results."
function radiplot(dO2_i, dtCO2_i, pfoc_i, psoc_i, proc_i)
    depths, dO2, dtCO2, pfoc, psoc, proc = radiwrap(dO2_i, dtCO2_i, pfoc_i,
        psoc_i, proc_i)
    ntps = size(dO2)[2]
    cmap = colormap("RdBu", ntps)
    cs = ntps
    p1 = plot(depths*100, dO2[:, 1]*1e3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="dO2")
    p2 = plot(depths*100, dtCO2[:, 1]*1e3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="dTCO2")
    p3 = plot(depths*100, pfoc[:, 1], legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="pF-POC")
    p4 = plot(depths*100, psoc[:, 1], legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="pS-POC")
    p5 = plot(depths*100, proc[:, 1], legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="pR-POC")
    for sp in 2:ntps
        plot!(p1, depths*100, dO2[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p2, depths*100, dtCO2[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p3, depths*100, pfoc[:, sp], legend=false, c=cmap[cs-sp+1])
        plot!(p4, depths*100, psoc[:, sp], legend=false, c=cmap[cs-sp+1])
        plot!(p5, depths*100, proc[:, sp], legend=false, c=cmap[cs-sp+1])
    end # for sp
    plot(p1, p2, p3, p3, p4, p5, layout=(3, 2))
    return depths, dO2, dtCO2, pfoc, psoc, proc
end # function radiplot

depths, dO2, dtCO2, pfoc, psoc, proc = radiplot(dO2_i, dtCO2_i,
    pfoc_i, psoc_i, proc_i)

end # module tst
