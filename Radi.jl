module Radi

using Colors, Plots
Plots.default(show=true)
Plots.closeall()
include("gsw_rho.jl")
include("Model.jl")

# Define model timesteps (all times in years)
stoptime = 50.0/128000  # how long to run for
interval = 5/128000  # duration of each model timestep
saveperXsteps = 1#28000  # save results in intervals of this many timesteps

# Import site-specific settings
include("IC_W29.jl")

"Convenient wrapper function for running RADI."
function wrap(
    dO2_i,
    dtCO2_i,
    dtNO3_i,
    dtSO4_i,
    dtPO4_i,
    dtNH4_i,
    dtH2S_i,
    dFeII_i,
    dMnII_i,
    pfoc_i,
    psoc_i,
    proc_i,
    pFeOH3_i,
    pMnO2_i,
)
    @time (
        depths,
        dO2,
        dtCO2,
        dtNO3,
        dtSO4,
        dtPO4,
        dtNH4,
        dtH2S,
        dFeII,
        dMnII,
        pfoc,
        psoc,
        proc,
        pFeOH3,
        pMnO2,
    ) = Model.timeloop(
        stoptime,
        interval,
        saveperXsteps,
        z_max,
        z_res,
        dbl,
        phiInf,
        phi0,
        beta,
        lambda_b,
        lambda_f,
        lambda_s,
        lambda_i,
        T,
        S,
        P,
        dO2_w,
        dtCO2_w,
        dtNO3_w,
        dtSO4_w,
        dtPO4_w,
        dtNH4_w,
        dtH2S_w,
        dFeII_w,
        dMnII_w,
        Fpom,
        Fpom_r,
        Fpom_s,
        Fpom_f,
        FMnO2,
        FFeOH3,
        rho_p,
        dO2_i,
        dtCO2_i,
        dtNO3_i,
        dtSO4_i,
        dtPO4_i,
        dtNH4_i,
        dtH2S_i,
        dFeII_i,
        dMnII_i,
        pfoc_i,
        psoc_i,
        proc_i,
        pFeOH3_i,
        pMnO2_i,
    )
end  # function radiwrap

"Run RADI and plot the results."
function profiles(
    dO2_i,
    dtCO2_i,
    dtNO3_i,
    dtSO4_i,
    dtPO4_i,
    dtNH4_i,
    dtH2S_i,
    dFeII_i,
    dMnII_i,
    pfoc_i,
    psoc_i,
    proc_i,
    pFeOH3_i,
    pMnO2_i,
)
    (
        depths,
        dO2,
        dtCO2,
        dtNO3,
        dtSO4,
        dtPO4,
        dtNH4,
        dtH2S,
        dFeII,
        dMnII,
        pfoc,
        psoc,
        proc,
        pFeOH3,
        pMnO2,
    ) = wrap(
        dO2_i,
        dtCO2_i,
        dtNO3_i,
        dtSO4_i,
        dtPO4_i,
        dtNH4_i,
        dtH2S_i,
        dFeII_i,
        dMnII_i,
        pfoc_i,
        psoc_i,
        proc_i,
        pFeOH3_i,
        pMnO2_i,
    )
    ntps = size(dO2)[2]
    cmap = colormap("RdBu", ntps)
    cs = ntps
    p1 = plot(depths*100, dO2[:, 1]*1e3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="O2")
    p2 = plot(depths*100, dtCO2[:, 1]*1e3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="TCO2")
    p8 = plot(depths*100, dtNO3[:, 1]*1e3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="dtNO3")
    p9 = plot(depths*100, dtSO4[:, 1]*1e3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="dtSO4")
    p10 = plot(depths*100, dtPO4[:, 1]*1e3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="dtPO4")
    p11 = plot(depths*100, dtNH4[:, 1]*1e3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="dtNH4")
    p12 = plot(depths*100, dtH2S[:, 1]*1e3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="dtH2S")
    p13 = plot(depths*100, dFeII[:, 1]*1e3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="dFeII")
    p14 = plot(depths*100, dMnII[:, 1]*1e3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="dMnII")
    p3 = plot(depths*100, pfoc[:, 1]*1e-3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="fPOC")
    p4 = plot(depths*100, psoc[:, 1]*1e-3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="sPOC")
    p5 = plot(depths*100, proc[:, 1]*1e-3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="rPOC")
    p6 = plot(depths*100, pFeOH3[:, 1]*1e-3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="pFeOH3")
    p7 = plot(depths*100, pMnO2[:, 1]*1e-3, legend=false, c=cmap[cs],
        xlabel="Depth / cm", ylabel="pMnO2")
    for sp in 2:ntps
        plot!(p1, depths*100, dO2[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p2, depths*100, dtCO2[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p8, depths*100, dtNO3[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p9, depths*100, dtSO4[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p10, depths*100, dtPO4[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p11, depths*100, dtNH4[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p12, depths*100, dtH2S[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p13, depths*100, dFeII[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p14, depths*100, dMnII[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p3, depths*100, pfoc[:, sp]*1e-3, legend=false, c=cmap[cs-sp+1])
        plot!(p4, depths*100, psoc[:, sp]*1e-3, legend=false, c=cmap[cs-sp+1])
        plot!(p5, depths*100, proc[:, sp]*1e-3, legend=false, c=cmap[cs-sp+1])
        plot!(p6, depths*100, pFeOH3[:, sp]*1e-3, legend=false, c=cmap[cs-sp+1])
        plot!(p7, depths*100, pMnO2[:, sp]*1e-3, legend=false, c=cmap[cs-sp+1])
    end  # for sp
    p0 = plot()
    plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p0, p0,
        layout=(4, 4))
    return (
        depths,
        dO2,
        dtCO2,
        dtNO3,
        dtSO4,
        dtPO4,
        dtNH4,
        dtH2S,
        dFeII,
        dMnII,
        pfoc,
        psoc,
        proc,
        pFeOH3,
        pMnO2,
    )
end  # function radiplot

Model.sayhello()
(
    depths,
    dO2,
    dtCO2,
    dtNO3,
    dtSO4,
    dtPO4,
    dtNH4,
    dtH2S,
    dFeII,
    dMnII,
    pfoc,
    psoc,
    proc,
    pFeOH3,
    pMnO2,
) = profiles(
    dO2_i,
    dtCO2_i,
    dtNO3_i,
    dtSO4_i,
    dtPO4_i,
    dtNH4_i,
    dtH2S_i,
    dFeII_i,
    dMnII_i,
    pfoc_i,
    psoc_i,
    proc_i,
    pFeOH3_i,
    pMnO2_i,
)

end  # module Radi
