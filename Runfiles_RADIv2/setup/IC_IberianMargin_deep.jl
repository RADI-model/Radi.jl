### Data taken from Epping et al., 2002 ###

# Define model depth steps (all depths in metres)
depthSed = 0.2  # depth of the sediment
z_res = 2e-3 #z stepsize
tspan = (0.0, 2000.0)
z_range = 0.0:z_res:depthSed
z = array_of_floats = collect(z_range)

# +
# wave
wave_height =1.0 #[m]
wave_period=6.1 #[s] typical wave period 
wavelength=(9.81*wave_period)/(2*pi) #[m] typical wavelength

depth=5.0 #[m] average depth of enclosed bay
# -

permeability=1e-13 # [m^2] typical sediment permeability for fine sediments

# Currents
U=0.02 #m/s

# Define sediment porosity parameters
phiInf = 0.74  # estimated sediment porosity at infinite depth from W-2
phi0 = 0.85  # estimated sediment porosity at the surface from W-2
beta = 33.0  # sediment porosity-depth relationship parameter

# Define characteristic depths
lambda_b = 0.08  # for bioturbation / m
lambda_f = 0.03  # for fast-degrading POC / m
lambda_s = 1.0  # for slow-degrading POC / m
lambda_i = 0.018  # for irrigation / m

# Define overlying water column properties
T = 4.354545  # (Epping et al., 2002)
S = 38.458  # practical salinity from GLODAP nearest station, bottom waters
P = 4380.0  # pressure at seafloor / dbar
rho_sw = gsw_rho(S, T, P) # seawater density / kg/m^3
# Concentrations all in mol/kg
dO2_w = 239.1e-6  # average dissolved oxygen (Epping et al., 2002)
dtCO2_w = 2313.8e-6 # total dissolved inorganic carbon from GLODAP nearest station, bottom waters
dtNO3_w = 19.636364e-6 # Epping et al., 2002
dtSO4_w = (29264.2e-6S/35)  # estimated computed from salinity (Millero, 2013)
dtPO4_w = 2.2165396e-6  # total phosphate from WOA nearest station, bottom waters
dtNH4_w = 0.255455e-6  # Epping et al., 2002
dtH2S_w = 0.0  #assumed
dCH4_w = 0.1e-6 #assumed
dFeII_w = 0.5e-9  # typical for deep-sea oxic bottom waters (Abadie et al., 2019)
dMnII_w = 0.5e-9  # typical for deep-sea oxic bottom waters (Morton et al., 2019)
dalk_w = 2582.2e-6  # total alkalinity from GLODAP nearest station, bottom waters
dSi_w = 9.66e-6  # total silicate from GLODAP nearest station, bottom waters

# Define organic matter flux to the surface sediment
Fpom = 10.0  # flux of POM to seafloor / g/m^2/a (Epping et al., 2002)
Fpom_r = 0.0  # refractory fraction of POM
Fpom_s = 0.36 # slow-degrading fraction of POM
Fpom_f = 0.64 # fast-degrading fraction of POM
FMnO2 = 0.0005  # typical for deep-sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
FFeOH3 = 0.0005  # typical for deep-sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
Fcalcite = 0.22  # flux of calcite to the seafloor / mol/m^2/a
Faragonite = 0.0  # flux of aragonite to the seafloor / mol/m^2/a
Fclay = (2.0 / 360.31)*5# flux of clay (montmorillonite) to the seafloor / mol/m^2/ to match sedimemtation rate of (Epping et al., 2002)
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


