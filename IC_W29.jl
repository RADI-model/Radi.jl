# Define model depth steps (all depths in metres)
z_res = 0.5e-2  # height of each depth step
z_max = 20e-2  # total height of the modelled sediment column

# Define diffusive boundary layer thickness / m
dbl = 1e-3

# Define sediment porosity parameters
phiInf = 0.74  # sediment porosity at infinite depth
phi0 = 0.85  # sediment porosity at the surface
beta = 33.0  # sediment porosity-depth relationship parameter

# Define characteristic depths
lambda_b = 0.08  # for bioturbation / m
lambda_f = 0.03  # for fast-degrading POC / m
lambda_s = 1.0  # for slow-degrading POC / m
lambda_i = 0.05  # for irrigation / m

# Define overlying water column properties
T = 1.4  # temperature / degC
S = 34.69  # practical salinity
P = 1.0  # pressure / dbar
rho_sw = gsw_rho(S, T, P) # seawater density / kg/m^3
dO2_w = 159.7e-6rho_sw  # dissolved oxygen / mol/m^3
dtCO2_w = 2324e-6rho_sw  # dissolved inorganic carbon / mol/m^3
dtPO4_w = 2.39e-6rho_sw  # total phosphate / mol/m^3

# Define organic matter flux to the surface sediment
Fpom = 36.45  # flux of POM to seafloor / g/m^2/a
Fpom_r = 0.15  # refractory fraction of POM
Fpom_s = 0.15  # slow-degrading fraction of POM
Fpom_f = 0.7  # fast-degrading fraction of POM
rho_pom = 2.65e6  # solid POM density / g/m^3

# Define initial conditions within the sediment (scalars or arrays)
dO2_i = dO2_w*2/3  # dissolved oxygen / mol/m^3
dtCO2_i = dtCO2_w*0.999  # dissolved inorganic carbon / mol/m^3
pfoc_i = 3e2  # fast-degrading particulate organic carbon / unit?
psoc_i = 3e3  # slow-degrading particulate organic carbon / unit?
proc_i = 3e4  # refractory particulate organic carbon / unit?
