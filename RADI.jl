module RADI

using Base.SimdLoop

include("gsw_rho.jl")

struct Solute
    now::Array{Float64,1}
    then::Array{Float64,1}
    above::Float64
    dvar::Array{Float64,1}
    save::Array{Float64,2}
end # struct Solute

struct Solid
    now::Array{Float64,1}
    then::Array{Float64,1}
    dvar::Array{Float64,1}
    save::Array{Float64,2}
end # struct Solute

SoluteOrSolid = SolidOrSolute = Union{Solid,Solute}

function Solute(var_start::Array{Float64,1}, above::Float64,
        dvar::Array{Float64,1}, var_save::Array{Float64,2})
    return Solute(copy(var_start), copy(var_start), above, dvar, var_save)
end # function Solute

function Solid(var_start::Array{Float64,1}, dvar::Array{Float64,1},
        var_save::Array{Float64,2})
    return Solid(copy(var_start), copy(var_start), dvar, var_save)
end # function Solid

"""
    timeprep(t_start, t_end, t_interval, saveperXsteps)

Prepare vectors of model timesteps and savepoints. All time units are in days.
"""
function timeprep(stoptime::Float64, interval::Float64, saveperXsteps::Int)
    timesteps::Array{Float64,1} = collect(0.0:interval:stoptime)
    ntps::Int = length(timesteps)
    savepoints::Array{Float64,1} = collect(1:saveperXsteps:ntps)
    # Save final timepoint if it's not already in the list
    if !(ntps in savepoints)
        append!(savepoints, ntps)
    end
    nsps::Int = length(savepoints)
    return timesteps, savepoints, ntps, nsps
end # function timeprep

"""
    model()

Run the RADI model.
"""
function model(stoptime::Float64, interval::Float64, saveperXsteps::Int,
    dO2_i, poc_i)
# ==============================================================================
# === User inputs/settings =====================================================
# ==============================================================================

# Model time grid: function inputs
# Model depth grid
z_res::Float64 = 2e-2 # 0.05e-2 # m
z_max::Float64 = 20e-2 # m

# Overlying water conditions
T::Float64 = 1.4 # temperature / degC
S::Float64 = 34.69 # practical salinity
rho_sw::Float64 = gsw_rho(S, T, 1) # seawater density [kg/m^3]
dO2_w::Float64 = 159.7e-6*rho_sw # dissolved oxygen [mol/m3]
po4_w::Float64 = 2.39e-6*rho_sw # phosphate from GLODAP at station location,
                                # bottom waters [mol/m3]
dbl::Float64 = 1e-3 # thickness at location taken from Sulpis et al. 2018 PNAS [m]

# Organic matter flux to the surface sediment
Ftot::Float64 = 36.45 # flux of particulate organic matter to seafloor [g/m2/a]
# Foc = 1.0 # flux of total organic carbon to the bottom [mol/m2/a]
rho_pom::Float64 = 2.65e6 # solid density [g/m3]

# Sediment porosity
# Porosity profile (porewater bulk fraction) fitted from station7 mooring3
# of cruise NBP98-2 by Sayles et al. DSR 2001
phiInf::Float64 = 0.74
phi0::Float64 = 0.85
beta::Float64 = 33.0

# Characteristic depths
lambda_b::Float64 = 0.08 # [m] characteristic depth of Archer et al. (2002),
                         # value here following Martin & Sayles (1990)
lambda_i::Float64 = 0.05 # [m] characteristic depth for irrigation

# ==============================================================================
# === The model: user changes nothing below here ===============================
# ==============================================================================

# Set up model time grid
timesteps::Array{Float64,1}, savepoints::Array{Float64,1}, ntps::Int,
    nsps::Int = timeprep(stoptime, interval, saveperXsteps)
sp::Int = 1

# Set up depth grid
z_res2::Float64 = z_res^2
depths::Array{Float64,1} = collect(-z_res:z_res:(z_max+z_res)) # in m
depths[1] = NaN
depths[end] = NaN
ndepths::Int = length(depths)

# "Redfield" ratios
RC::Float64 = @. 1.0/(6.9e-3po4_w/(1e-6rho_sw) + 6e-3)
# ^P:C computed as a function of SRP from Galbraith and Martiny PNAS 2015
RN::Float64 = 11.0 # value at 60 degS from Martiny et al. Nat G 2013
RP::Float64 = 1.0 # Redfield ratio for P in the deep sea
M_CH2O::Float64 = 30.031 # g/mol
M_NH3::Float64 = 17.031 # g/mol
M_H3PO4::Float64 = 97.994 # g/mol
# M_OM = @. M_CH2O + (RN/RC)*M_NH3 + (RP/RC)*M_H3PO4
# # ^ Organic Matter molar mass [g of OM per mol of OC]
M_POM::Float64 = RC*M_CH2O + RN*M_NH3 + RP*M_H3PO4 # molar mass of POM in g/mol
Ftot_mol::Float64 = Ftot/M_POM
Foc::Float64 = Ftot_mol*RC

# Depth-dependent porosity
phi::Array{Float64,1} = @. (phi0 - phiInf)*exp(-beta*depths) + phiInf
phiS::Array{Float64,1} = 1.0 .- phi # solid volume fraction
phiS_phi::Array{Float64,1} = phiS./phi
tort2::Array{Float64,1} = @. 1.0 - 2.0log(phi)
# ^tortuosity squared from Boudreau (1996, GCA)
delta_phi::Array{Float64,1} = @. -beta*(phi0 - phiInf)*exp(-beta*depths)
# delta_phi[2] = 0.0 # don't do this
delta_phiS::Array{Float64,1} = -delta_phi
delta_tort2i::Array{Float64,1} = @. 2.0delta_phi/(phi*tort2^2)
delta_tort2_tort2::Array{Float64,1} = delta_tort2i.*tort2

# Bioturbation (for solids)
D_bio_0::Float64 = @. 0.0232e-4*(1e2Foc)^0.85
# ^[m2/a] surf bioturb coeff, Archer et al. (2002)
D_bio::Array{Float64,1} = @. D_bio_0*exp(-(depths/lambda_b)^2)*
    dO2_w/(dO2_w + 0.02)
# ^[m2/a] bioturb coeff, Archer et al (2002)
delta_D_bio::Array{Float64,1} = @. -2.0depths*D_bio/lambda_b^2

# Organic matter degradation parameters
krefractory::Array{Float64,1} = @. 80.25D_bio_0*exp(-depths)
# ^[/a] from Archer et al (2002)

# Solid fluxes and solid initial conditions
# Ftot = Foc.*M_OM # total sediment flux [g/m2/a]
x0::Float64 = Ftot/(rho_pom*phiS[2])
# ^[m/a] bulk burial velocity at sediment-water interface
xinf::Float64 = x0*phiS[2]/phiS[end-1]
# ^[m/a] bulk burial velocity at the infinite depth
u::Array{Float64,1} = xinf*phi[end-1]./phi # [m/a] porewater burial velocity
w::Array{Float64,1} = xinf*phiS[end-1]./phiS # [m/a] solid burial velocity

# Biodiffusion depth-attenuation: see Boudreau (1996); Fiadeiro and Veronis
# (1977)
Peh::Array{Float64,1} = @. w*z_res/2.0D_bio
# ^one half the cell Peclet number (Eq. 97 in Boudreau 1996)
# when Peh<<1, biodiffusion dominates, when Peh>>1, advection dominates
sigma::Array{Float64,1} = @. 1.0/tanh(Peh) - 1.0/(Peh) # Eq. 96 in Boudreau 1996
# sigma[2] = 0.0 # don't do this
sigma1m::Array{Float64,1} = 1.0 .- sigma
sigma1p::Array{Float64,1} = 1.0 .+ sigma

# vvv NOT YET IN THE PARAMETERS PART OF THE DOCUMENTATION vvvvvvvvvvvvvvvvvvvvvv
# Temperature-dependent "free solution" diffusion coefficients
D_dO2::Float64 = @. 0.034862 + 0.001409T
# ^[m2/a] oxygen diffusion coefficient from Li and Gregory
D_dO2_tort2::Array{Float64,1} = D_dO2./tort2

# Irrigation (for solutes)
alpha_0::Float64 = @. 11.0*(atan((1e2Foc*5.0 - 400.0)/400.0)/pi + 0.5) - 0.9 +
    20.0*(dO2_w/(dO2_w + 0.01))*exp(-dO2_w/0.01)*1e2Foc/(1e2Foc + 30.0)
# ^[/a] from Archer et al (2002)
alpha::Array{Float64,1} = @. alpha_0*exp(-(depths/lambda_i)^2)
# ^[/a] Archer et al (2002) the depth of 5 cm was changed

APPW::Array{Float64,1} = @. w - delta_D_bio - delta_phiS*D_bio/phiS
TR::Float64 = 2.0z_res*tort2[2]/dbl
Foc_phiS_0::Float64 = Foc/phiS[2]
zr_Db_0::Float64 = 2.0z_res/D_bio[2]
# ^^^ NOT YET IN THE PARAMETERS PART OF THE DOCUMENTATION ^^^^^^^^^^^^^^^^^^^^^^

"Prepare Solute with a constant start value."
function makeSolute(var_start::Float64, above::Float64, D_var::Array{Float64})
    var_start = fill(var_start, ndepths)
    var_start[1] = NaN
    var_start[end] = NaN
    var_save = fill(NaN, (ndepths-2, nsps+1))
    var_save[:, 1] = var_start[2:end-1]
    return Solute(var_start, above, D_var, var_save)
end # function makeSolute

"Prepare Solute with starting array provided."
function makeSolute(var_start::Array{Float64,1}, above::Float64,
        D_var::Array{Float64})
    var_start = vcat(NaN, var_start, NaN)
    var_save = fill(NaN, (ndepths-2, nsps+1))
    var_save[:, 1] = var_start[2:end-1]
    return Solute(var_start, above, D_var, var_save)
end # function makeSolute

"Prepare Solid with a constant start value."
function makeSolid(var_start::Float64, D_var::Array{Float64})
    var_start = fill(var_start, ndepths)
    var_start[1] = NaN
    var_start[end] = NaN
    var_save = fill(NaN, (ndepths-2, nsps+1))
    var_save[:, 1] = var_start[2:end-1]
    return Solid(var_start, D_var, var_save)
end # function makeSolid

"Prepare Solid with starting array provided."
function makeSolid(var_start::Array{Float64,1}, D_var::Array{Float64})
    var_start = vcat(NaN, var_start, NaN)
    var_save = fill(NaN, (ndepths-2, nsps+1))
    var_save[:, 1] = var_start[2:end-1]
    return Solid(var_start, D_var, var_save)
end # function makeSolid

# Create variables to model
dO2 = makeSolute(dO2_i, dO2_w, D_dO2_tort2)
poc = makeSolid(poc_i, D_bio)

"Calculate the above-surface value for a solute."
function surfacesolute(then::Array{Float64,1}, above::Float64)
    # # Equation following Boudreau (1996, method-of-lines):
    # n = 2 # ambiguous value from Eq. (104)
    # return then[3] + (above - then[2])*2z_res/(dbl*phi[2]^(n+1))
    # Or, equation following RADI-Matlab and CANDI-Fortran:
    return then[3] + (above - then[2])*TR
end # function surfacesolute

"Calculate the above-surface value for a solid."
function surfacesolid(then::Array{Float64,1})
    return then[3] + (Foc_phiS_0 - w[2]*then[2])*zr_Db_0
end # function surfacesolid

"Calculate the below-bottom value for a solid or solute."
function bottom(then::Array{Float64,1})
    return then[end-2]
end # function bottom

"Substitute in the above-surface and below-bottom values for a Solute."
function substitute!(var::Solute)
    var.then[1] = surfacesolute(var.then, var.above)
    var.then[end] = bottom(var.then)
end # function substitute!

"Substitute in the above-surface and below-bottom values for a Solid."
function substitute!(var::Solid)
    var.then[1] = surfacesolid(var.then)
    var.then[end] = bottom(var.then)
end # function substitute!

"React a Solute or Solid."
function react!(z::Int, var::SolidOrSolute, rate::Float64)
    var.now[z] += interval*rate
end # function react!

"Calculate advection rate for a solute."
function advectsolute(then_z1p::Float64, then_z1m::Float64, u_z::Float64,
        delta_phi_z::Float64, phi_z::Float64, delta_tort2_tort2_z::Float64,
        D_var::Float64)
    return -(u_z - delta_phi_z*D_var/phi_z -
        D_var*delta_tort2_tort2_z)*(then_z1p - then_z1m)/(2.0z_res)
end # function advect

"Advect a Solute."
function advect!(z::Int, var::Solute)
    var.now[z] += interval*advectsolute(var.then[z+1], var.then[z-1], u[z],
        delta_phi[z], phi[z], delta_tort2_tort2[z], var.dvar[z])
end # function advect!

"Calculate advection rate for a solid."
function advectsolid(then_z::Float64, then_z1p::Float64, then_z1m::Float64,
        APPW_z::Float64, sigma_z::Float64, sigma1p_z::Float64,
        sigma1m_z::Float64)
    return -APPW_z*(sigma1m_z*then_z1p + 2.0sigma_z*then_z -
        sigma1p_z*then_z1m)/(2.0z_res)
end # function advectsolid

"Advect a Solid."
function advect!(z::Int, var::Solid)
    # var[z] += -interval*APPW[z]*(sigma1m[z]*var0[z+1] + 2.0sigma[z]*var0[z] -
    #     sigma1p[z]*var0[z-1])/(2.0z_res)
    var.now[z] += interval*advectsolid(var.then[z], var.then[z+1],
        var.then[z-1], APPW[z], sigma[z], sigma1p[z], sigma1m[z])
end # function advect!

"Calculate diffusion rate of a solute or solid."
function diffuse(then_z1m::Float64, then_z::Float64, then_z1p::Float64,
        D_var::Float64)
    return (then_z1m - 2.0then_z + then_z1p)*D_var/z_res2
end # function diffuse

"Diffuse a Solute or Solid."
function diffuse!(z::Int, var::SolidOrSolute)
    var.now[z] += interval*diffuse(var.then[z-1], var.then[z], var.then[z+1],
        var.dvar[z])
end # function diffuse!

"Calculate irrigation rate of a solute."
function irrigate(then_z::Float64, var_w::Float64, alpha_z::Float64)
    return alpha_z*(var_w - then_z)
end # function irrigate

"Irrigate a Solute throughout the sediment."
function irrigate!(z::Int, var::Solute)
    var.now[z] += interval*irrigate(var.then[z], var.above, alpha[z])
end # function irrigate!

# ===== Run RADI run: main model loop ==========================================
for t in 1:ntps
    tsave::Bool = t in savepoints # i.e. do we save this time?
    # Substitutions above and below the modelled sediment column
    substitute!(dO2)
    substitute!(poc)
    @simd for z in 2:(ndepths-1)
    # ~~~ BEGIN SEDIMENT PROCESSING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # --- First, do all the physical processes -----------------------------
        # Dissolved oxygen (solute)
        advect!(z, dO2)
        diffuse!(z, dO2)
        irrigate!(z, dO2)
        # Particulate organic carbon (solid)
        advect!(z, poc)
        diffuse!(z, poc)
        # --- Now do the reactions! --------------------------------------------
        # Calculate maximum reaction rates based on previous timestep
        R_poc::Float64 = -poc.then[z]*krefractory[z]
        R_dO2::Float64 = R_poc*phiS_phi[z]
        # Check maximum reaction rates are possible after other processes have
        # acted in this timestep, and correct them if not
        if dO2.now[z] + interval*R_dO2 < 0.0
            R_dO2 = -dO2.now[z]/interval
            R_poc = R_dO2/phiS_phi[z]
        end # if
        if poc.now[z] + interval*R_poc < 0.0
            R_poc = -poc.now[z]/interval
            R_dO2 = R_poc*phiS_phi[z]
        end # if
        react!(z, dO2, R_dO2)
        react!(z, poc, R_poc)
    # ~~~ END SEDIMENT PROCESSING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Save output if we are at a savepoint
        if tsave
            dO2.save[z-1, sp+1] = dO2.now[z]
            poc.save[z-1, sp+1] = poc.now[z]
            if z == ndepths-1
                println("RADI: reached savepoint $sp (step $t of $ntps)...")
                sp += 1
            end # if
        end # if
    end # for z in 2:(ndepths-1)
    # Copy results into "previous step" arrays
    @simd for z in 2:(ndepths-1)
        dO2.then[z] = dO2.now[z]
        poc.then[z] = poc.now[z]
    end # for z in 2:(ndepths-1)
end # for t
# ===== End of main model loop =================================================
return depths[2:end-1], dO2.save, poc.save
end # function model

say_RADI() = println("RADI done!")

end # module RADI
