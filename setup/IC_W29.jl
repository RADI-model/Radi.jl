# Set model run name
const modelrun = "IC_W29"

# Define model timesteps (all times in years)
const stoptime = 10.0  # how long to run for
const interval = 1/128000  # duration of each model timestep
const saveperXsteps = 128000  # save results in intervals of this many timesteps

# Define model depth steps (all depths in metres)
const z_res = 0.2e-2  # height of each depth step
const z_max = 20e-2  # total height of the modelled sediment column

# Define diffusive boundary layer thickness / m
const dbl = 1e-3

# Define sediment porosity parameters
const phiInf = 0.74  # sediment porosity at infinite depth
const phi0 = 0.85  # sediment porosity at the surface
const beta = 33.0  # sediment porosity-depth relationship parameter

# Define characteristic depths
const lambda_b = 0.08  # for bioturbation / m
const lambda_f = 0.03  # for fast-degrading POC / m
const lambda_s = 1.0  # for slow-degrading POC / m
const lambda_i = 0.05  # for irrigation / m

# Define overlying water column properties
const T = 1.4  # temperature / degC
const S = 34.69  # practical salinity
const P = 4380.0  # pressure at seafloor / dbar
const rho_sw = gsw_rho(S, T, P) # seawater density / kg/m^3
# Concentrations all in mol/m^3
const dO2_w = 159.7e-6rho_sw  # dissolved oxygen
const dtCO2_w = 2324e-6rho_sw  # dissolved inorganic carbon
const dtNO3_w = 36.93e-6rho_sw  # nitrate from GLODAP at station location, bottom waters
const dtSO4_w = (29264.2e-6S/35)rho_sw  # estimated omputed from salinity (Millero, 2013)
const dtPO4_w = 2.39e-6rho_sw  # total phosphate
const dtNH4_w = 0.0  # assumed
const dtH2S_w = 0.0  # assumed
const dFeII_w = 0.5e-9rho_sw  # typical for deep-sea oxic bottom waters (Abadie et al., 2019)
const dMnII_w = 0.5e-9rho_sw  # typical for deep-sea oxic bottom waters (Morton et al., 2019)
const dalk_w = 2426e-6rho_sw  # total alkalinity from GLODAP at station location, bottom waters
const dCa_w = 0.02128 / 40.087 * S / 1.80655 * rho_sw  # calcium from salinity (RT67 via PyCO2SYS)
const dSi_w = 120e-6rho_sw  # total silicate

# Define organic matter flux to the surface sediment
const Fpom = 6.561  # flux of POM to seafloor / g/m^2/a
const Fpom_r = 0.03  # refractory fraction of POM
const Fpom_s = 0.27  # slow-degrading fraction of POM
const Fpom_f = 0.7  # fast-degrading fraction of POM
const FMnO2 = 0.0005  # typical for deep-sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
const FFeOH3 = 0.0005  # typical for deep-sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
const Fcalcite = 0.22  # flux of calcite to the seafloor / mol/m^2/a
const Faragonite = 0.0  # flux of aragonite to the seafloor / mol/m^2/a
const Fclay = 2.0 / 360.31  # flux of clay (montmorillonite) to the seafloor / mol/m^2/a
const rho_p = 2.65e6  # average density of all solid matter / g/m^3

# Define initial conditions within the sediment (scalars or arrays)
const dO2_i = copy(dO2_w)  # dissolved oxygen / mol/m^3
const dtCO2_i = copy(dtCO2_w)  # dissolved inorganic carbon / mol/m^3
const dtNO3_i = copy(dtNO3_w)
const dtSO4_i = copy(dtSO4_w)
const dtPO4_i = copy(dtPO4_w)
const dtNH4_i = copy(dtNH4_w)
const dtH2S_i = copy(dtH2S_w)
const dFeII_i = copy(dFeII_w)
const dMnII_i = copy(dMnII_w)
const dalk_i = copy(dalk_w)
const dCa_i = copy(dCa_w)
const pfoc_i = 0.0  # fast-degrading particulate organic carbon / unit?
const psoc_i = 3e2  # slow-degrading particulate organic carbon / unit?
const proc_i = 6e2  # refractory particulate organic carbon / unit?
const pFeOH3_i = 0.0
const pMnO2_i = 0.0
const pcalcite_i = 4.5e3
const paragonite_i = 0.0
const pclay_i = 1.0e3
