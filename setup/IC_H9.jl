# Set model run name
const modelrun = "IC_H9"

# Define model depth steps (all depths in metres)
depthSed = 0.20  # depth of the sediment
z_res = 5e-3 #z stepsize
tspan = (0.0, 20.0) # in years
z_range = 0.0:z_res:depthSed
z = array_of_floats = collect(z_range)

# wave (not currently included)
wave_height =1.0 #[m]
wave_period=6.1 #[s] typical wave period 
wavelength=(9.81*wave_period)/(2*pi) #[m] typical wavelength

depth=5.0 #[m] average depth of enclosed bay

# bottom water currents
U= 0.02 #m/s

# Define sediment parameters
phiInf = 0.74 # sediment porosity at infinite depth
phi0 = 0.91  # sediment porosity at the surface
beta = 33.0  # sediment porosity-depth relationship parameter

permeability=6.6e-12 # sediment permeability

# Define characteristic depths
lambda_b = 0.08  # for bioturbation / m
lambda_f = 0.03  # for fast-degrading POC / m
lambda_s = 1.0  # for slow-degrading POC / m
lambda_i = 0.018 # for irrigation / m

# Define overlying water column properties
T = 2.2  # temperature / degC
S = 34.9  # practical salinity
P = 5312.4  # pressure at seafloor / dbar

# Concentrations all in mol/kg (will be made mol/m3 in the model for ensemble runs where S and T might differ) 
dO2_w = 266.6e-6  # dissolved oxygen from GLODAP at station location, bottom waters
dtCO2_w = 2186e-6  # dissolved inorganic carbon from GLODAP at station location, bottom waters
dtNO3_w = 20.0668e-6  # nitrate from GLODAP at station location, bottom waters
dtSO4_w = (29180e-6S/35)  # estimated omputed from salinity (Millero, 2013)
dtPO4_w = 1.3561e-6  # total phosphate from GLODAP at station location, bottom waters
dtNH4_w = 0.0  # assumed
dtH2S_w = 0.0  # assumed
dFeII_w = 0.5e-9  # typical for deep-sea oxic bottom waters (Abadie et al., 2019)
dMnII_w = 0.5e-9  # typical for deep-sea oxic bottom waters (Morton et al., 2019)
dalk_w = 2342e-6  # total alkalinity from GLODAP at station location, bottom waters
dSi_w = 120e-6  # total silicate

# Define organic matter flux to the surface sediment
Fpom = 6.0  # flux of POM to seafloor / g/m^2/a
Fpom_r = 0.03  # refractory fraction of POM
Fpom_s = 0.27  # slow-degrading fraction of POM
Fpom_f = 0.70  # fast-degrading fraction of POM
FMnO2 = 0.0005  # typical for deep-sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
FFeOH3 = 0.0005  # typical for deep-sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
Fcalcite = 0.2  # flux of calcite to the seafloor / mol/m^2/a
Faragonite = 0.0  # flux of aragonite to the seafloor / mol/m^2/a
Fclay = 26.0 / 360.31  # flux of clay (montmorillonite) to the seafloor / mol/m^2/a
rho_p = 2.65e6  # average density of all solid matter / g/m^3

# Define initial conditions within the sediment (scalars or arrays)
dO2_i = copy(dO2_w)  # dissolved oxygen
dtCO2_i = copy(dtCO2_w)  # dissolved inorganic carbon
dtNO3_i = copy(dtNO3_w)
dtSO4_i = copy(dtSO4_w)
dtPO4_i = copy(dtPO4_w)
dtNH4_i = copy(dtNH4_w)
dtH2S_i = copy(dtH2S_w)
dFeII_i = copy(dFeII_w)
dMnII_i = copy(dMnII_w)
dalk_i = copy(dalk_w)
dCa_i = copy(dCa_w)
pfoc_i = 0.0  # fast-degrading particulate organic carbon / unit?
psoc_i = 3.0e2  # slow-degrading particulate organic carbon / unit?
proc_i = 6.0e2  # refractory particulate organic carbon / unit?
pFeOH3_i = 0.0
pMnO2_i = 0.0
pcalcite_i = 4.5e3
paragonite_i = 0.0
pclay_i = 1.0e3
