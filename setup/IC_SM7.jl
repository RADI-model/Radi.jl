# Set model run name
const modelrun = "IC_SM7"

# Define model timesteps (all times in years)
const stoptime = 10.0  # how long to run for
const interval = 1/128000  # duration of each model timestep
const saveperXsteps = 128000  # save results in intervals of this many timesteps

# Define model depth steps (all depths in metres)
const z_res = 0.2e-2  # height of each depth step
const z_max = 40e-2  # total height of the modelled sediment column

# Define diffusive boundary layer thickness / m
const dbl = 0.715e-3

# Define sediment porosity parameters
const phiInf = 0.87  # sediment porosity at infinite depth
const phi0 = 0.91  # sediment porosity at the surface
const beta = 33.0  # sediment porosity-depth relationship parameter

# Define characteristic depths
const lambda_b = 0.08  # for bioturbation / m
const lambda_f = 0.03  # for fast-degrading POC / m
const lambda_s = 1.0  # for slow-degrading POC / m
const lambda_i = 0.05  # for irrigation / m

# Define overlying water column properties
const T = 0.84  # temperature / degC
const S = 34.696  # practical salinity
const P = 3932.8  # pressure at seafloor / dbar
const rho_sw = gsw_rho(S, T, P) # seawater density / kg/m^3
# Concentrations all in mol/m^3
const dO2_w = 215.7e-6rho_sw  # dissolved oxygen from GLODAP at station location, bottom waters
const dtCO2_w = 2260e-6rho_sw  # dissolved inorganic carbon from GLODAP at station location, bottom waters
const dtNO3_w = 32.2416e-6rho_sw  # nitrate from GLODAP at station location, bottom waters
const dtSO4_w = (29264.2e-6S/35)rho_sw  # estimated computed from salinity (Millero, 2013)
const dtPO4_w = 2.2428e-6rho_sw  # total phosphate from GLODAP at station location, bottom waters
const dtNH4_w = 0.0  # assumed
const dtH2S_w = 0.0  # assumed
const dFeII_w = 0.5e-9rho_sw  # typical for deep-sea oxic bottom waters (Abadie et al., 2019)
const dMnII_w = 0.5e-9rho_sw  # typical for deep-sea oxic bottom waters (Morton et al., 2019)
const dalk_w = 2365e-6rho_sw  # total alkalinity from GLODAP at station location, bottom waters
const dCa_w = 0.02128 / 40.087 * S / 1.80655 * rho_sw  # calcium from salinity (RT67 via PyCO2SYS)
const dSi_w = 120e-6rho_sw  # total silicate

# Define organic matter flux to the surface sediment
const Fpom = 4.6242672045283015  # flux of POM to seafloor / g/m^2/a
const Fpom_r = 0.03  # refractory fraction of POM
const Fpom_s = 0.27  # slow-degrading fraction of POM
const Fpom_f = 0.70  # fast-degrading fraction of POM
const FMnO2 = 0.0005  # typical for deep-sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
const FFeOH3 = 0.0005  # typical for deep-sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
const Fcalcite = 0.25  # flux of calcite to the seafloor / mol/m^2/a
const Faragonite = 0.0  # flux of aragonite to the seafloor / mol/m^2/a
const Fclay = 32.0 / 360.31  # flux of clay (montmorillonite) to the seafloor / mol/m^2/a
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
const psoc_i = 3.0e2  # slow-degrading particulate organic carbon / unit?
const proc_i = 6.0e2  # refractory particulate organic carbon / unit?
const pFeOH3_i = 0.0
const pMnO2_i = 0.0
const pcalcite_i = 4.5e3
const paragonite_i = 0.0
const pclay_i = 1.0e3
