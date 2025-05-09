### Data taken from Epping et al., 2002 ###

# Define model depth steps (all depths in metres)
depthSed = 0.2  # depth of the sediment
z_res = 0.2e-2 #z stepsize
tspan = (0.0, 4000.0)
z_range = 0.0:z_res:depthSed
z = array_of_floats = collect(z_range)

# +
# wave
wave_height =1.0 #[m]
wave_period=6.1 #[s] typical wave period 
wavelength=(9.81*wave_period)/(2*pi) #[m] typical wavelength

depth=5.0 #[m] average depth of enclosed bay
# -

permeability=5e-12 # [m^2] typical sediment permeability for sandy sediments

# Currents
U=0.02 #m/s

# Define sediment porosity parameters
phiInf = 0.45 # estimated sediment porosity based on the data from the North Sea
phi0 = 0.71  # estimated sediment porosity based on the data from the North Sea
beta = 33.0  # sediment porosity-depth relationship parameter

# Define characteristic depths
lambda_b = 0.08  # for bioturbation / m
lambda_f = 0.03  # for fast-degrading POC / m
lambda_s = 1.0  # for slow-degrading POC / m
lambda_i = 0.05  # for irrigation / m

# Define overlying water column properties
T = 12.108333  # (Epping et al., 2002)
S = 38.458  # practical salinity from GLODAP nearest station, bottom waters
P = 4380.0  # pressure at seafloor / dbar
rho_sw = gsw_rho(S, T, P) # seawater density / kg/m^3
# Concentrations all in mol/kg
dO2_w = 216.333333e-6  # average dissolved oxygen (Epping et al., 2002)
dtCO2_w = 2313.8e-6 # total dissolved inorganic carbon from GLODAP nearest station, bottom waters
dtNO3_w = 10.800000e-6 # Epping et al., 2002
dtSO4_w = (29264.2e-6S/35)  # estimated computed from salinity (Millero, 2013)
dtPO4_w = 2.2165396e-6  # total phosphate from WOA nearest station, bottom waters
dtNH4_w = 0.325833e-6  # Epping et al., 2002
dtH2S_w = 0.0  #assumed
dCH4_w = 0.1e-6 #assumed
dFeII_w =  4.186537551421602e-8 # mean bottom water concentration North Sea (Siems et al., 2024)
dMnII_w = 2.1788197604572425e-7 # mean bottom water concentration North Sea (Siems et al., 2024)
dalk_w = 2582.2e-6  # total alkalinity from GLODAP nearest station, bottom waters
dSi_w = 9.66e-6  # total silicate from GLODAP nearest station, bottom waters

# Define organic matter flux to the surface sediment
Fpom = 44.215395617325186  # flux of POM to seafloor / g/m^2/a (Epping et al., 2002)
Fpom_r = 0.0  # refractory fraction of POM
Fpom_s = 0.27348635644808994# slow-degrading fraction of POM
Fpom_f = 0.7265136435519101# fast-degrading fraction of POM
FMnO2 = 0.02 * (365.25/1000)  # North Sea flux from (Slomp et al., 1996) converted to mol/m^2/a
FFeOH3 = 0.10 * (365.25/1000)# North Sea flux from (Slomp et al., 1996) converted to mol/m^2/a
Fcalcite = 0.22  # flux of calcite to the seafloor / mol/m^2/a
Faragonite = 0.0  # flux of aragonite to the seafloor / mol/m^2/a
Fclay = (2.0 / 360.31)*300# flux of clay (montmorillonite) to the seafloor / mol/m^2/ to match sedimemtation rate of (Epping et al., 2002)
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


