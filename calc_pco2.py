import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import PyCO2SYS as pyco2

mpl.rcParams["font.family"] = "Open Sans"
mpl.rcParams["font.size"] = mpl.rcParams["axes.titlesize"] = 10

# Set general conditions
npts = 100
pH_guess = np.linspace(7.0, 9.0, npts)
alk = np.full(npts, 2300.0) * 1e-6
dic = np.full(npts, 2101.0) * 1e-6
tempC = np.full(npts, 15.0)
pres = np.full(npts, 0.0)
psal = np.full(npts, 35.0)
silicate = np.full(npts, 0.0) * 1e-6
phosphate = np.full(npts, 0.0) * 1e-6
ammonia = np.full(npts, 0.0) * 1e-6
sulfide = np.full(npts, 0.0) * 1e-6
whichKs = np.full(npts, 10)
whoseTB = np.full(npts, 2)
whoseKSO4 = np.full(npts, 1)
whoseKF = np.full(npts, 1)
whichR = np.full(npts, 1)
pHscale = np.full(npts, 1)

# Assemble intermediates
totals = pyco2.salts.assemble(
    psal, silicate, phosphate, ammonia, sulfide, whichKs, whoseTB
)
ks = pyco2.equilibria.assemble(
    tempC, pres, totals, pHscale, whichKs, whoseKSO4, whoseKF, whichR
)
free2tot = pyco2.convert.free2tot(totals["TSO4"], ks["KSO4"])

# Get carbonate alkalinity
alk_parts = pyco2.solve.get.AlkParts(dic, pH_guess, free2tot, totals, ks)
alk_carb = (
    alk
    - alk_parts["BAlk"]
    - alk_parts["OH"]
    + alk_parts["Hfree"]
    - alk_parts["PAlk"]
    - alk_parts["SiAlk"]
)

# Solve for pH, calc_pco2 style
gamma = dic / alk_carb
b_h = ks["K1"] * (1 - gamma)
discr = b_h ** 2 - 4 * ks["K1"] * ks["K2"] * (1 - 2 * gamma)
h_new = (np.sqrt(discr) - b_h) / 2
pH_new = -np.log10(h_new)

# Get true pH
pH_true = pyco2.solve.get.pHfromTATC(alk, dic, totals, ks)

# Visualise
fig, ax = plt.subplots(dpi=300)
ax.plot(pH_guess, pH_guess, c="xkcd:ocean blue", label="Initial guess")
ax.plot(pH_guess, pH_new, c="xkcd:reddish orange", label="Final estimate")
ax.plot(pH_guess, pH_true, c="k", label="True value")
ax.set_aspect("equal")
ax.set_xlim([7, 9])
ax.set_ylim([7, 9])
ax.set_xlabel("Initial pH guess")
ax.set_ylabel("pH")
ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
ax.grid(alpha=0.3)
ax.legend(edgecolor="k")
plt.tight_layout()
plt.savefig("figures/calc_pco2.png")
