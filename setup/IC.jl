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
U= 0.07 #m/s

# Define sediment parameters
phiInf = 0.45 # sediment porosity at infinite depth
phi0 = 0.71  # sediment porosity at the surface
beta = 33.0  # sediment porosity-depth relationship parameter

permeability=6.6e-12 # sediment permeability

# Define characteristic depths
lambda_b = 0.08  # for bioturbation / m
lambda_f = 0.03  # for fast-degrading POC / m
lambda_s = 1.0  # for slow-degrading POC / m
lambda_i = 0.018 # for irrigation / m

# Define overlying water column properties
T = 7.967995901806943  # temperature / degC
S = 35.136  # practical salinity
P = 167.97593908354008  # pressure at seafloor / dbar (1atm = 10.1325 dbar )

# Concentrations all in mol/kg (will be made mol/m3 in the model for ensemble runs where S and T might differ) 
dO2_w = 250e-6 # dissolved oxygen
dtCO2_w = 2164.0e-6 # dissolved inorganic carbon
dtNO3_w = 4.4e-6 # dissolved nitrate
dtSO4_w = (29264.2e-6S/35) # estimated omputed from salinity (Millero, 2013)
dtPO4_w = 1.0e-6 # total phosphate
dtNH4_w = 6.0e-6
dtH2S_w = 0.0  
dFeII_w =  4.186537551421602e-8 
dMnII_w = 2.1788197604572425e-7 
dCH4_w = 0.1e-6
dalk_w = 2300.4e-6  # total alkalinity
dSi_w = 4.1e-6 # total silicate

# Define organic matter flux to the surface sediment
Fpom = 135.0 # flux of POM to seafloor / g/m^2/a, divide by Mpom in Model to get mol/m^2/a
Fpom_r = 0.00  # refractory fraction of POM
Fpom_s = 0.48068700143424825 # slow-degrading fraction of POM
Fpom_f = 0.5193129985657517 # fast-degrading fraction of POM
FMnO2 = 0.02 * (365.25/1000)  # flux of manganese oxide to the seafloor / mol/m^2/a
FFeOH3 = 0.10 * (365.25/1000) # flux of iron oxide to the seafloor / mol/m^2/a
Fcalcite = 0.5  # flux of calcite to the seafloor / mol/m^2/a
Faragonite = 0.0  # flux of aragonite to the seafloor / mol/m^2/a
Fclay = (2.0 / 360.31)*80 # currently unreactive, can be used to play with the sedimentation rate
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
dCH4_i = copy(dCH4_w)
dalk_i = copy(dalk_w)
pfoc_i = 0.0  # fast-degrading particulate organic carbon / unit?
psoc_i = 200  # slow-degrading particulate organic carbon / unit?
proc_i = 6e2  # refractory particulate organic carbon / unit?
pFeOH3_i = 0.0
pMnO2_i = 0.0
pcalcite_i = 4.5e3
paragonite_i = 0.0
pclay_i = 1.0e3
