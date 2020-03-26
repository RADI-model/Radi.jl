module RADI

include("gsw_rho.jl")

"""
    timeprep(t_start, t_end, t_interval, saveperXsteps)

Prepare vectors of model timesteps and savepoints. All time units are in days.
"""
function timeprep(stoptime, interval, saveperXsteps)
    timesteps = 0.0:interval:stoptime
    ntps = length(timesteps)
    savepoints = collect(1:saveperXsteps:ntps)
    # Save final timepoint if it's not already in the list
    if !(ntps in savepoints)
        append!(savepoints, ntps)
    end
    nsps = length(savepoints)
    return timesteps, savepoints, ntps, nsps
end # function timeprep

"""
    model()

Run the RADI model.
"""
function model(stoptime, interval, saveperXsteps)
# ==============================================================================
# === User inputs/settings =====================================================
# ==============================================================================
# Overlying water conditions
T = 1.4 # temperature / degC
S = 34.69 # practical salinity
rho_sw = gsw_rho(S, T, 1) # seawater density [kg/m^3]
oxy_w = 159.7e-6*rho_sw # dissolved oxygen [mol/m3]
po4_w = 2.39e-6*rho_sw # phosphate from GLODAP at station location, bottom
                       # waters [mol/m3]
# Model depth grid
z_res = 0.001
z_max = 0.2
# Sediment porosity
# Porosity profile (porewater bulk fraction) fitted from station7 mooring3
# of cruise NBP98-2 by Sayles et al. DSR 2001
phi0 = 0.85
phiInf = 0.74
beta = 33.0
# Organic matter flux to the surface sediment
Foc = 1.0 # flux of total organic carbon to the bottom [mol/m2/a]
rho_pom = 2.65e6 # solid density [g/m3]
poc_initial = 100.0 # initial POC in the sediment
# Diffusive boundary layer
dbl = 1e-3 # thickness at location taken from Sulpis et al. 2018 PNAS [m]
# Characteristic depths
lambda_b = 0.08 # [m] characteristic depth of Archer et al. (2002), value
                # here following Martin & Sayles (1990)
lambda_i = 0.05 # [m] characteristic depth for irrigation

# ==============================================================================
# === The model: user changes nothing below here ===============================
# ==============================================================================

# Set up model time variables
timesteps, savepoints, ntps, nsps = timeprep(stoptime, interval, saveperXsteps)
sp = 1

# "Redfield" ratios
RC = @. 1/(6.9e-3po4_w/(1e-6rho_sw) + 6e-3)
# ^P:C computed as a function of SRP from Galbraith and Martiny PNAS 2015
RN = 11 # value at 60 degS from Martiny et al. Nat G 2013
RP = 1 # Redfield ratio for P in the deep sea
M_OM = @. 30.03 + (RN/RC)*17.03 + (RP/RC)*97.994
# ^ Organic Matter molar mass [g of OM per mol of OC]

# Create arrays for modelled variables
z_res2 = z_res^2
depths = collect(-z_res:z_res:(z_max+z_res)) # in m
depths[1] = NaN
depths[end] = NaN
ndepths = length(depths)
oxy0 = fill(oxy_w, ndepths)
oxy = copy(oxy0)
poc0 = fill(poc_initial, ndepths)
poc = copy(poc0)

# Depth-dependent porosity
phi = @. (phi0 - phiInf)*exp(-beta*depths) + phiInf
phiS = 1.0 .- phi # solid volume fraction
phiS_phi = phiS./phi
tort2 = @. 1.0 - 2.0log(phi) # tortuosity squared from Boudreau (1996, GCA)
tort = sqrt.(tort2) # tortuosity

# Solid fluxes and solid initial conditions
Ftot = Foc.*M_OM # total sediment flux [g/m2/a]
v0 = Ftot/(rho_pom*phiS[2]) # [m/a] bulk burial velocity at sediment-water
                            # interface
vinf = v0*phiS[2]/phiS[end-1] # [m/a] bulk burial velocity at the infinite depth
u = vinf*phi[end-1]./phi # [m/a] porewater burial velocity
w = vinf*phiS[end-1]./phiS # [m/a] solid burial velocity

# Temperature-dependent "free solution" diffusion coefficients
D_oxy = @. 0.034862 + 0.001409T # [m2/a] oxygen diffusion coefficient from Li
                                # and Gregory
D_oxy_tort2 = D_oxy./tort2

# Bioturbation (for solids)
D_bio_0 = @. 0.0232e-4*(1e2Foc)^0.85 # [m2/a] surf bioturb coeff, Archer et al.
                                     # (2002)
D_bio = @. D_bio_0*exp(-(depths/lambda_b)^2)*oxy_w/(oxy_w + 0.02)
# ^[m2/a] bioturb coeff, Archer et al (2002)

# Irrigation (for solutes)
alpha_0 = @. 11.0*(atan((1e2Foc*5.0 - 400.0)/400.0)/pi + 0.5) - 0.9 +
    20.0*(oxy_w/(oxy_w + 0.01))*exp(-oxy_w/0.01)*1e2Foc/(1e2Foc + 30.0)
# ^[/a] from Archer et al (2002)
alpha = @. alpha_0*exp(-(depths/lambda_i)^2)
# ^[/a] Archer et al (2002) the depth of 5 cm was changed

# Depth-dependent porosity and diffusion coefficient loss
delta_phi = @. -beta*(phi0 - phiInf)*exp(-beta*depths)
delta_phiS = @. beta*(phi0 - phiInf)*exp(-beta*depths)
delta_tort2 = @. 2.0beta*(phi0 - phiInf)/
    (phi0 + phiInf*(exp(beta*depths) - 1.0))
delta_D_bio = @. -(2.0D_bio_0*depths/lambda_b^2)*exp(-(depths/lambda_b)^2)*
    oxy_w/(oxy_w + 0.02)
# delta_phi = append!([0.0], diff(phi)) # depth-dependent porosity loss
# delta_phiS = append!([0.0], diff(phiS)) # depth-dependent solid fraction gain
# delta_tort2 = append!([0.0], diff(tort.^2)) # depth-dependent tortuosity gain
# delta_D_bio = append!([0.0], diff(D_bio)) # [m/a]

# Biodiffusion depth-attenuation: see Boudreau (1996); Fiadeiro and Veronis
# (1977)
Peh = @. w*z_res/2.0D_bio # one half the cell Peclet number (Eq. 97 in Boudreau
                          # 1996)
# when Peh<<1, biodiffusion dominates, when Peh>>1, advection dominates
sigma = @. 1.0/tanh(Peh) - 1.0/(Peh) # Eq. 96 in Boudreau 1996

# Organic matter degradation parameters
krefractory = @. 80.25D_bio_0*exp(-depths) # [/a] from Archer et al (2002)

# Short-cut transport variables
APPW = @. w - delta_D_bio - delta_phiS*D_bio/phiS
DFF = @. (tort2*delta_phi/phi - delta_tort2)/tort2^2
# DBF = @. D_bio*(phiS - delta_phi) # not used?
# TR = @. 2z_res*tort2/dbl

# Create arrays for saved variables
oxy_save = fill(NaN, (ndepths-2, nsps))
poc_save = fill(NaN, (ndepths-2, nsps))

"Substitute in above-surface value for solutes."
function surface!(var0::Array{Float64,1}, var_w::Float64)
    # # Equation following Boudreau (1996, method-of-lines):
    # n = 2 # ambiguous value from Eq. (104)
    # var0[1] = var0[3] + (var_w - var0[2])*2z_res/(dbl*phi[2]^(n+1))
    # Or, equation following RADI-Matlab and CANDI-Fortran:
    var0[1] = var0[3] + (var_w - var0[2])*2.0z_res*tort2[2]/dbl
end # function surface!

"Substitute in above-surface value for solids."
function surface!(var0::Array{Float64,1}, Fvar::Float64, D_var_0::Float64)
    var0[1] = var0[3] + (Fvar/phiS[2] - w[2]*var0[2])*2.0z_res/D_var_0
end # function surface!

"Substitute in below-bottom value."
function bottom!(var0::Array{Float64,1})
    var0[end] = var0[end-2]
end # function bottom!

"Substitute in above-surface and below-bottom values for solutes."
function substitute!(var0::Array{Float64,1}, var_w::Float64)
    surface!(var0, var_w)
    bottom!(var0)
end # function substitute!

"Substitute in above-surface and below-bottom values for solids."
function substitute!(var0::Array{Float64,1}, Fvar::Float64,
        D_var_0::Float64)
    surface!(var0, Fvar, D_var_0)
    bottom!(var0)
end # function substitute!

"Reactions throughout the sediment, solutes and solids."
function react(rate::Float64)
    return interval*rate
end # function react!

"Reactions throughout the sediment, solutes and solids."
function react!(z::Int, var::Array{Float64,1}, rate::Float64)
    var[z] += react(rate)
end # function react!



"Diffusion throughout the sediment, solutes and solidst."
function diffuse!(z::Int, var0::Array{Float64,1}, var::Array{Float64,1},
        D_var::Float64)
    var[z] += interval*(var0[z-1] - 2.0var0[z] + var0[z+1])*D_var/z_res2
end # function diffuse!

"Irrigation of solutes only throughout the sediment."
function irrigate!(z::Int, var0::Array{Float64,1}, var::Array{Float64,1},
        var_w::Float64)
    var[z] += interval*alpha[z]*(var_w - var0[z])
end # function irrigate!

# ===== Run RADI run: main model loop ==========================================
for t in 1:ntps
    tsave = t in savepoints # i.e. do we save this time?

    # Substitutions above and below the modelled sediment column
    substitute!(oxy0, oxy_w)
    substitute!(poc0, Foc, D_bio_0)

    for z in 2:(ndepths-1)
    # ~~~ BEGIN SEDIMENT PROCESSING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # First, do all the physical processes
        diffuse!(z, oxy0, oxy, D_oxy_tort2[z])
        irrigate!(z, oxy0, oxy, oxy_w)

        diffuse!(z, poc0, poc, D_bio[z])

        # Calculate maximum reaction rates based on previous timestep
        R_poc = -poc0[z]*krefractory[z]
        R_oxy = R_poc*phiS_phi[z]

        # Check maximum reaction rates are possible after other processes
        # have acted in this timestep, and correct them if not
        if oxy[z] + react(R_oxy) < 0.0
            R_oxy = -oxy[z]/interval
            R_poc = R_oxy/phiS_phi[z]
        end # if
        if poc[z] + react(R_poc) < 0.0
            R_poc = -poc[z]/interval
            R_oxy = R_poc*phiS_phi[z]
        end # if

        # Now do the reactions!
        react!(z, oxy, R_oxy)
        react!(z, poc, R_poc)
    # ~~~ END SEDIMENT PROCESSING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Save output if we are at a savepoint
        if tsave
            if t == 1
                oxy_save[z-1, sp] = oxy0[z]
                poc_save[z-1, sp] = poc0[z]
            else
                oxy_save[z-1, sp] = oxy[z]
                poc_save[z-1, sp] = poc[z]
            end
            if z == ndepths-1
                println("RADI: reached savepoint $sp (step $t of $ntps)...")
                sp += 1
            end # if
        end # if
    end # for z in 2:(ndepths-1)

    # Copy results into "previous step" arrays
    for z in 2:(ndepths-1)
        oxy0[z] = oxy[z]
        poc0[z] = poc[z]
    end # for z in 2:(ndepths-1)
end # for t
# ===== End of main model loop =================================================
return depths[2:end-1], oxy_save, poc_save
end # function model

say_RADI() = println("RADI done!")

end # module RADI
