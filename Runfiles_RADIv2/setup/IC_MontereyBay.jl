# Set model run name
const modelrun = "IC_MontereyBay.jl"

### Data taken from Berelson et al., 2003 ###

# Define model depth steps (all depths in metres)
depthSed = 0.2  # depth of the sediment
z_res = 2e-3 #z stepsize
tspan = (0.0, 1500.0)
z_range = 0.0:z_res:depthSed
z = array_of_floats = collect(z_range)

# +
# wave
wave_height =1.0 #[m]
wave_period=6.1 #[s] typical wave period 
wavelength=(9.81*wave_period)/(2*pi) #[m] typical wavelength

depth=5.0 #[m] average depth of enclosed bay
# -

permeability=6.6e-12 # [m^2] typical sediment permeability for sandy sediments
# Currents
U=0.02 #m/s

# Define sediment porosity parameters
phiInf = 0.50  # sediment porosity at infinite depth
phi0 = 0.70  # sediment porosity at the surface
beta = 33.0  # sediment porosity-depth relationship parameter

# Define characteristic depths
lambda_b = 0.08  # for bioturbation / m
lambda_f = 0.03  # for fast-degrading POC / m
lambda_s = 1.0  # for slow-degrading POC / m
lambda_i = 0.05  # for irrigation / m

# Define overlying water column properties
T = 6.0  # between 4-8 temperature / degC
S = 35.136  # practical salinity
P = 100.0  # pressure at seafloor / dbar
# Concentrations all in mol/kg
dO2_w = 154.42857142857142e-6  # average dissolved oxygen (Berelson et al., 2003)
dtCO2_w = 2160.3e-6 # total dissolved inorganic carbon from GLODAP nearest station, bottom waters
dtNO3_w = 22.5e-6 # Berelson et al., 2003
dtSO4_w = (29264.2e-6S/35)  # estimated omputed from salinity (Millero, 2013)
dtPO4_w = 2.39e-6  # total phosphate (Berelson et al., 2003)
dtNH4_w = 3.0e-6  # assumed
dtH2S_w = 0.0  #assumed
dCH4_w = 0.1e-6 #assumed
dFeII_w = 0.5e-9  # typical for deep-sea oxic bottom waters (Abadie et al., 2019)
dMnII_w = 0.5e-9  # typical for deep-sea oxic bottom waters (Morton et al., 2019)
dalk_w = 2230.7e-6  # total alkalinity from GLODAP nearest station, bottom waters
dSi_w = 30e-6  # total silicate (Berelson et al., 2003)

# Define organic matter flux to the surface sediment
Fpom = 122.2  # flux of POM to seafloor / g/m^2/a
Fpom_r = 0.00  # refractory fraction of POM
Fpom_s = 0.48391178424171155 # slow-degrading fraction of POM
Fpom_f = 0.5160882157582884 # fast-degrading fraction of POM
FMnO2 = 0.02 * (365.25/1000)  # North Sea self flux from (Slomp et al., 1996) converted to mol/m^2/a
FFeOH3 = 0.10 * (365.25/1000)# North Sea shelf flux from (Slomp et al., 1996) converted to mol/m^2/a
Fcalcite = 0.22  # flux of calcite to the seafloor / mol/m^2/a
Faragonite = 0.0  # flux of aragonite to the seafloor / mol/m^2/a
Fclay = (2.0 / 360.31)*350 # flux of clay (montmorillonite) to the seafloor / mol/m^2/a
rho_p = 2.65e6  # average density of all solid matter / g/m^3

# Define initial conditions within the sediment (scalars or arrays)
dO2_i = copy(dO2_w)  # dissolved oxygen / mol/m^3
dtCO2_i = copy(dtCO2_w)  # dissolved inorganic carbon / mol/m^3
dtNO3_i = copy(dtNO3_w)
dtSO4_i = copy(dtSO4_w)
dtPO4_i = copy(dtPO4_w)
dtNH4_i = copy(dtNH4_w)
dtH2S_i = copy(dtH2S_w)
dFeII_i = copy(dFeII_w)
dMnII_i = copy(dMnII_w)
dCH4_i = copy(dCH4_w)
dalk_i = copy(dalk_w)
pfoc_i = 0.0  # fast-degrading particulate organic carbon / unit?
psoc_i = 3e2  # slow-degrading particulate organic carbon / unit?
proc_i = 6e2  # refractory particulate organic carbon / unit?
pFeOH3_i = 0.0
pMnO2_i = 0.0
pcalcite_i = 4.5e3
paragonite_i = 0.0
pclay_i = 1.0e3


