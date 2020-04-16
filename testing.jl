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

# Define model depth steps (all depths in metres)
z_res = 0.5e-2 # height of each depth step
z_max = 20e-2 # total height of the modelled sediment column

# Define diffusive boundary layer thickness in metres
dbl = 1e-3

# Define sediment porosity parameters
phiInf = 0.74 # sediment porosity at infinite depth
phi0 = 0.85 # sediment porosity at the surface
beta = 33.0 # sediment porosity-depth relationship parameter

# Define characteristic depths
lambda_b = 0.08 # for bioturbation
lambda_i = 0.05 # for irrigation

# Define overlying water column properties
T = 1.4 # temperature / degC
S = 34.69 # practical salinity
P = 1.0 # pressure / dbar
rho_sw = gsw_rho(S, T, P) # seawater density / kg/m^3
dO2_w = 159.7e-6*rho_sw # dissolved oxygen / mol/m^3
dtPO4_w = 2.39e-6*rho_sw # total phosphate / mol/m^3

# Define organic matter flux to the surface sediment
Fpom = 36.45 # flux of POM to seafloor / g/m^2/a
rho_pom = 2.65e6 # solid POM density / g/m^3

# Define initial conditions within the sediment (scalars or arrays)
dO2_i = dO2_w*2/3 # dissolved oxygen / mol/m^3
poc_i = 0.0 # particulate organic carbon / unit?

# Collect model arguments except input conditions
RADIargs = (stoptime, interval, saveperXsteps, z_max, z_res, dbl, phiInf, phi0,
    beta, lambda_b, lambda_i, T, S, P, dO2_w, dtPO4_w, Fpom, rho_pom)

function radiplot(dO2_i, poc_i)
    @time depths, oxyx, pocx = RADI.model(RADIargs..., dO2_i, poc_i)
    oxy = oxyx.save
    poc = pocx.save
    ntps = size(oxy)[2]
    cmap = colormap("RdBu", ntps)
    cs = ntps
    p1 = plot(depths*100, oxy[:, 1]*1e3, legend=false, c=cmap[cs])
    p2 = plot(depths*100, poc[:, 1], legend=false, c=cmap[cs])
    for sp in 2:ntps
        plot!(p1, depths*100, oxy[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p2, depths*100, poc[:, sp], legend=false, c=cmap[cs-sp+1])
    end # for sp
    plot(p1, p2, layout=(2, 1))
    return depths, oxy, poc
end # function radiplot

showprofile = false
if showprofile
    Profile.clear()
    @profile depths, oxy, poc = RADI.model(RADIargs..., dO2_i, poc_i)
    ProfileView.view()
    RADI.say_RADI()
else
    depths, oxy, poc = radiplot(dO2_i, poc_i)
end

end # module tst
