module Model

using Base.SimdLoop
include("gsw_rho.jl")
include("Params.jl")
include("React.jl")

"Define Solute type."
struct Solute
    now::Array{Float64,1}
    then::Array{Float64,1}
    above::Float64
    dvar::Array{Float64,1}
    save::Array{Float64,2}
end  # struct Solute

"Define Solid type."
struct Solid
    now::Array{Float64,1}
    then::Array{Float64,1}
    above::Float64
    dvar::Array{Float64,1}
    save::Array{Float64,2}
end  # struct Solute

"Constructor for a Solute."
function Solute(var_start::Array{Float64,1}, above::Float64,
        dvar::Array{Float64,1}, var_save::Array{Float64,2})
    return Solute(copy(var_start), copy(var_start), above, dvar, var_save)
end  # function Solute

"Constructor for a Solid."
function Solid(var_start::Array{Float64,1}, above::Float64,
        dvar::Array{Float64,1}, var_save::Array{Float64,2})
    return Solid(copy(var_start), copy(var_start), above, dvar, var_save)
end  # function Solid

SoluteOrSolid = SolidOrSolute = Union{Solid,Solute}
FloatOrArray = ArrayOrFloat = Union{Float64,Array{Float64,1}}

"Prepare vectors of model timesteps and savepoints. All time units are in days."
function preptime(stoptime::Float64, interval::Float64, saveperXsteps::Int)
    timesteps::Array{Float64,1} = collect(0.0:interval:stoptime)
    ntps::Int = length(timesteps)
    savepoints::Array{Int64,1} = collect((1+saveperXsteps):saveperXsteps:ntps)
    # Save final timepoint too if it's not already in the list
    if !(ntps in savepoints)
        append!(savepoints, ntps)
    end
    nsps::Int = length(savepoints)
    return timesteps, savepoints, ntps, nsps
end  # function preptime

"Prepare the model's depth vector. All depth units are in metres."
function prepdepth(depth_res::Float64, depth_max::Float64)
    depths = collect(-depth_res:depth_res:(depth_max+depth_res))  # in m
    depths[1] = NaN
    depths[end] = NaN
    ndepths = length(depths)
    depth_res2 = depth_res^2
    return depths, ndepths, depth_res2
end  # function prepdepth

"Assemble depth-dependent porosity parameters."
function porosity(phi0::Float64, phiInf::Float64, beta::Float64,
        depths::Array{Float64})
    phi = Params.phi(phi0, phiInf, beta, depths)
    phiS = Params.phiS(phi)  # solid volume fraction
    phiS_phi = phiS./phi
    tort2 = Params.tort2(phi)  # tortuosity squared from Boudreau (1996, GCA)
    delta_phi = Params.delta_phi(phi0, phiInf, beta, depths)
    delta_phiS = Params.delta_phiS(delta_phi)
    delta_tort2i = Params.delta_tort2i(delta_phi, phi, tort2)
    delta_tort2i_tort2 = delta_tort2i.*tort2
    return phi, phiS, phiS_phi, tort2, delta_phi, delta_phiS, delta_tort2i_tort2
end  # function porosity

"Define 'Redfield' ratios and OM stoichiometry."
function stoichiometry(T::Float64, S::Float64, P::Float64, dtPO4_w::Float64,
        Fpom::Float64)
    rho_sw = gsw_rho(S, T, P)  # seawater density [kg/m^3]
    RC, RN, RP = Params.redfield(dtPO4_w, rho_sw)
    Mpom = Params.rmm_pom(RC, RN, RP)
    Fpom_mol = Fpom/Mpom
    Fpoc = Fpom_mol*RC
end  # function stoichiometry

"Run the iterative Radi model."
function timeloop(
    stoptime::Float64,
    interval::Float64,
    saveperXsteps::Int,
    z_max::Float64,
    z_res::Float64,
    dbl::Float64,
    phiInf::Float64,
    phi0::Float64,
    beta::Float64,
    lambda_b::Float64,
    lambda_f::Float64,
    lambda_s::Float64,
    lambda_i::Float64,
    T::Float64,
    S::Float64,
    P::Float64,
    dO2_w::Float64,
    dtCO2_w::Float64,
    dtNO3_w::Float64,
    dtSO4_w::Float64,
    dtPO4_w::Float64,
    dtNH4_w::Float64,
    dtH2S_w::Float64,
    dFeII_w::Float64,
    dMnII_w::Float64,
    Fpom::Float64,
    Fpom_r::Float64,
    Fpom_s::Float64,
    Fpom_f::Float64,
    FMnO2::Float64,
    FFeOH3::Float64,
    rho_p::Float64,
    dO2_i::FloatOrArray,
    dtCO2_i::FloatOrArray,
    dtNO3_i::FloatOrArray,
    dtSO4_i::FloatOrArray,
    dtPO4_i::FloatOrArray,
    dtNH4_i::FloatOrArray,
    dtH2S_i::FloatOrArray,
    dFeII_i::FloatOrArray,
    dMnII_i::FloatOrArray,
    pfoc_i::FloatOrArray,
    psoc_i::FloatOrArray,
    proc_i::FloatOrArray,
    pFeOH3_i::FloatOrArray,
    pMnO2_i::FloatOrArray,
)

println("Radi preparing to run...")

# Set up model time and depth grids
timesteps, savepoints, ntps, nsps = preptime(stoptime, interval, saveperXsteps)
depths, ndepths, z_res2 = prepdepth(z_res, z_max)
sp = 1  # initialise savepoints

# Calculate depth-dependent porosity
phi, phiS, phiS_phi, tort2, delta_phi, delta_phiS, delta_tort2i_tort2 =
    porosity(phi0, phiInf, beta, depths)

# Define 'Redfield' ratios and OM stoichiometry
rho_sw = gsw_rho(S, T, P)  # seawater density [kg/m^3]
RC, RN, RP = Params.redfield(dtPO4_w, rho_sw)
Mpom = Params.rmm_pom(RC, RN, RP)  # g/mol
Fpom_mol = Fpom/Mpom  # mol/m^2/a
Fpoc = Fpom_mol*RC  # mol/m^2/a
# Split total flux into fast-slow-refractory portions
Ffoc = Fpoc*Fpom_f
Fsoc = Fpoc*Fpom_s
Froc = Fpoc*Fpom_r
if Fpom_f + Fpom_s + Fpom_r != 1.0
    println("\nRadi WARNING: the fractions of POM do not add up to 1!\n")
end
# `Fp` = total sediment flux to bottom in g/m^2/a
M_MnO2 = 86.9368  # g/mol
M_FeOH3 = 106.867  # g/mol
Fp = Fpom + FMnO2*M_MnO2 + FFeOH3*M_FeOH3

# Bioturbation (for solids)
D_bio_0 = Params.D_bio_0(Fpoc)
# ^[m2/a] surf bioturb coeff, Archer et al. (2002)
D_bio = Params.D_bio(depths, D_bio_0, lambda_b, dO2_w)
# ^[m2/a] bioturb coeff, Archer et al (2002)
delta_D_bio = Params.delta_D_bio(depths, D_bio, lambda_b)

# Organic matter degradation parameters
krefractory = Params.krefractory(depths, D_bio_0)
kfast = Params.kfast(Fpoc, depths, lambda_f)
kslow = Params.kslow(Fpoc, depths, lambda_s)
# ^[/a] from Archer et al (2002)

# Solid fluxes and solid initial conditions
x0 = Params.x0(Fp, rho_p, phiS[2])
# ^[m/a] bulk burial velocity at sediment-water interface
xinf = Params.xinf(x0, phiS[2], phiS[end-1])
# ^[m/a] bulk burial velocity at the infinite depth
u = Params.u(xinf, phi)  # [m/a] porewater burial velocity
w = Params.w(xinf, phiS)  # [m/a] solid burial velocity

# Biodiffusion depth-attenuation: see Boudreau (1996); Fiadeiro & Veronis (1977)
Peh = Params.Peh(w, z_res, D_bio)
# ^one half the cell Peclet number (Eq. 97 in Boudreau 1996)
# When Peh<<1, biodiffusion dominates, when Peh>>1, advection dominates
sigma = Params.sigma(Peh)
sigma1m = 1.0 .- sigma
sigma1p = 1.0 .+ sigma

# vvv NOT YET IN THE PARAMETERS PART OF THE DOCUMENTATION vvvvvvvvvvvvvvvvvv
# Temperature-dependent "free solution" diffusion coefficients
D_dO2 = Params.D_dO2(T)
D_dtCO2 = Params.D_dtCO2(T)
D_dtNO3 = Params.D_dtNO3(T)
D_dtSO4 = Params.D_dtSO4(T)
D_dtPO4 = Params.D_dtPO4(T)
D_dtNH4 = Params.D_dtNH4(T)
D_dtH2S = Params.D_dtH2S(T)
D_dMnII = Params.D_dMn(T)
D_dFeII = Params.D_dFe(T)
D_dO2_tort2 = D_dO2./tort2
D_dtCO2_tort2 = D_dtCO2./tort2
D_dtNO3_tort2 = D_dtNO3./tort2
D_dtSO4_tort2 = D_dtSO4./tort2
D_dtPO4_tort2 = D_dtPO4./tort2
D_dtNH4_tort2 = D_dtNH4./tort2
D_dtH2S_tort2 = D_dtH2S./tort2
D_dMnII_tort2 = D_dMnII./tort2
D_dFeII_tort2 = D_dFeII./tort2

# Irrigation (for solutes)
alpha_0 = Params.alpha_0(Fpoc, dO2_w)  # [/a] from Archer et al (2002)
alpha = Params.alpha(alpha_0, depths, lambda_i)  # [/a] Archer et al (2002)

# Miscellaneous convenience variables
APPW = Params.APPW(w, delta_D_bio, delta_phiS, D_bio, phiS)
TR = Params.TR(z_res, tort2[2], dbl)
zr_Db_0 = 2.0z_res/D_bio[2]
# ^^^ NOT YET IN THE PARAMETERS PART OF THE DOCUMENTATION ^^^^^^^^^^^^^^^^^^

"Prepare Solute with a constant start value."
function makeSolute(var_start::Float64, above::Float64, D_var::Array{Float64})
    var_start = fill(var_start, ndepths)
    var_start[1] = NaN
    var_start[end] = NaN
    var_save = fill(NaN, (ndepths-2, nsps+1))
    var_save[:, 1] = var_start[2:end-1]
    return Solute(var_start, above, D_var, var_save)
end  # function makeSolute

"Prepare Solute with starting array provided."
function makeSolute(var_start::Array{Float64,1}, above::Float64,
        D_var::Array{Float64})
    var_start = vcat(NaN, var_start, NaN)
    var_save = fill(NaN, (ndepths-2, nsps+1))
    var_save[:, 1] = var_start[2:end-1]
    return Solute(var_start, above, D_var, var_save)
end  # function makeSolute

"Prepare Solid with a constant start value."
function makeSolid(var_start::Float64, above::Float64, D_var::Array{Float64})
    var_start = fill(var_start, ndepths)
    var_start[1] = NaN
    var_start[end] = NaN
    above_phiS_0 = above/phiS[2]
    var_save = fill(NaN, (ndepths-2, nsps+1))
    var_save[:, 1] = var_start[2:end-1]
    return Solid(var_start, above_phiS_0, D_var, var_save)
end  # function makeSolid

"Prepare Solid with starting array provided."
function makeSolid(var_start::Array{Float64,1}, above::Float64,
        D_var::Array{Float64})
    var_start = vcat(NaN, var_start, NaN)
    above_phiS_0 = above/phiS[2]
    var_save = fill(NaN, (ndepths-2, nsps+1))
    var_save[:, 1] = var_start[2:end-1]
    return Solid(var_start, above_phiS_0, D_var, var_save)
end  # function makeSolid

"Calculate the above-surface value for a solute."
function surfacesolute(then::Array{Float64,1}, above::Float64)
    # # Equation following Boudreau (1996, method-of-lines):
    # n = 2 # ambiguous value from Eq. (104)
    # return then[3] + (above - then[2])*2z_res/(dbl*phi[2]^(n+1))
    # Or, equation following Radi-Matlab and CANDI-Fortran:
    return then[3] + (above - then[2])*TR
end  # function surfacesolute

"Calculate the above-surface value for a solid."
function surfacesolid(then::Array{Float64,1}, above::Float64)
    return then[3] + (above - w[2]*then[2])*zr_Db_0
end  # function surfacesolid

"Calculate the below-bottom value for a solid or solute."
function bottom(then::Array{Float64,1})
    return then[end-2]
end  # function bottom

"Substitute in the above-surface and below-bottom values for a Solute."
function substitute!(var::Solute)
    var.then[1] = surfacesolute(var.then, var.above)
    var.then[end] = bottom(var.then)
end  # function substitute!

"Substitute in the above-surface and below-bottom values for a Solid."
function substitute!(var::Solid)
    var.then[1] = surfacesolid(var.then, var.above)
    var.then[end] = bottom(var.then)
end  # function substitute!

"React a Solute or Solid."
function react!(var::SolidOrSolute, z::Int, rate::Float64)
    change = interval * rate
    if change > var.now[z]
        var.now[z] = 0.0
        println("Warning: overriding negative concentration with zero.")
    else
        var.now[z] += change
    end  # if else
end  # function react!

"Calculate advection rate for a solute."
function advectsolute(then_z1p::Float64, then_z1m::Float64, u_z::Float64,
        delta_phi_z::Float64, phi_z::Float64, delta_tort2i_tort2_z::Float64,
        D_var::Float64)
    return -(u_z - delta_phi_z*D_var/phi_z -
        D_var*delta_tort2i_tort2_z)*(then_z1p - then_z1m)/(2.0z_res)
end  # function advect

"Advect a Solute."
function advect!(var::Solute, z::Int)
    var.now[z] += interval*advectsolute(var.then[z+1], var.then[z-1], u[z],
        delta_phi[z], phi[z], delta_tort2i_tort2[z], var.dvar[z])
end  # function advect!

"Calculate advection rate for a solid."
function advectsolid(then_z::Float64, then_z1p::Float64, then_z1m::Float64,
        APPW_z::Float64, sigma_z::Float64, sigma1p_z::Float64,
        sigma1m_z::Float64)
    return -APPW_z*(sigma1m_z*then_z1p + 2.0sigma_z*then_z -
        sigma1p_z*then_z1m)/(2.0z_res)
end  # function advectsolid

"Advect a Solid."
function advect!(var::Solid, z::Int)
    var.now[z] += interval*-APPW[z]*(sigma1m[z]*var.then[z+1] +
        2.0sigma[z]*var.then[z] - sigma1p[z]*var.then[z-1])/(2.0z_res)
# No idea why the approach below is so much slower, only for this function?!
    # var.now[z] += interval*advectsolid(var.then[z], var.then[z+1],
    #     var.then[z-1], APPW[z], sigma[z], sigma1p[z], sigma1m[z])
end  # function advect!

"Calculate diffusion rate of a solute or solid."
function diffuse(then_z1m::Float64, then_z::Float64, then_z1p::Float64,
        D_var::Float64)
    return (then_z1m - 2.0then_z + then_z1p)*D_var/z_res2
end  # function diffuse

"Diffuse a Solute or Solid."
function diffuse!(var::SolidOrSolute, z::Int)
    var.now[z] += interval*diffuse(var.then[z-1], var.then[z], var.then[z+1],
        var.dvar[z])
end  # function diffuse!

"Calculate irrigation rate of a solute."
function irrigate(then_z::Float64, above::Float64, alpha_z::Float64)
    return alpha_z*(above - then_z)
end  # function irrigate

"Irrigate a Solute throughout the sediment."
function irrigate!(var::Solute, z::Int)
    var.now[z] += interval*irrigate(var.then[z], var.above, alpha[z])
end  # function irrigate!

# ===== Run Radi run! ==========================================================
# Create variables to model
dO2 = makeSolute(dO2_i, dO2_w, D_dO2_tort2)
dtCO2 = makeSolute(dtCO2_i, dtCO2_w, D_dtCO2_tort2)
dtNO3 = makeSolute(dtNO3_i, dtNO3_w, D_dtNO3_tort2)
dtSO4 = makeSolute(dtSO4_i, dtSO4_w, D_dtSO4_tort2)
dtPO4 = makeSolute(dtPO4_i, dtPO4_w, D_dtPO4_tort2)
dtNH4 = makeSolute(dtNH4_i, dtNH4_w, D_dtNH4_tort2)
dtH2S = makeSolute(dtH2S_i, dtH2S_w, D_dtH2S_tort2)
dFeII = makeSolute(dFeII_i, dFeII_w, D_dFeII_tort2)
dMnII = makeSolute(dMnII_i, dMnII_w, D_dMnII_tort2)
pfoc = makeSolid(pfoc_i, Ffoc, D_bio)
psoc = makeSolid(psoc_i, Fsoc, D_bio)
proc = makeSolid(proc_i, Froc, D_bio)
pFeOH3 = makeSolid(pFeOH3_i, FFeOH3, D_bio)
pMnO2 = makeSolid(pMnO2_i, FMnO2, D_bio)
# Main Radi model loop
for t in 1:ntps
    tsave = t in savepoints  # i.e. do we save after this step?
    # Substitutions above and below the modelled sediment column
    substitute!(dO2)
    substitute!(dtCO2)
    substitute!(dtNO3)
    substitute!(dtSO4)
    substitute!(dtPO4)
    substitute!(dtNH4)
    substitute!(dtH2S)
    substitute!(dFeII)
    substitute!(dMnII)
    substitute!(pfoc)
    substitute!(psoc)
    substitute!(proc)
    substitute!(pFeOH3)
    substitute!(pMnO2)
    @simd for z in 2:(ndepths-1)
    # ~~~ BEGIN SEDIMENT PROCESSING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # --- First, do all the physical processes -----------------------------
        # Dissolved oxygen (solute)
        advect!(dO2, z)
        diffuse!(dO2, z)
        irrigate!(dO2, z)
        # Dissolved inorganic carbon (solute)
        advect!(dtCO2, z)
        diffuse!(dtCO2, z)
        irrigate!(dtCO2, z)
        # Nitrate (solute)
        advect!(dtNO3, z)
        diffuse!(dtNO3, z)
        irrigate!(dtNO3, z)
        # Sulfate (solute)
        advect!(dtSO4, z)
        diffuse!(dtSO4, z)
        irrigate!(dtSO4, z)
        # Phosphate (solute)
        advect!(dtPO4, z)
        diffuse!(dtPO4, z)
        irrigate!(dtPO4, z)
        # Ammonium (solute)
        advect!(dtNH4, z)
        diffuse!(dtNH4, z)
        irrigate!(dtNH4, z)
        # Hydrogen sulfide (solute)
        advect!(dtH2S, z)
        diffuse!(dtH2S, z)
        irrigate!(dtH2S, z)
        # Iron-II (solute)
        advect!(dFeII, z)
        diffuse!(dFeII, z)
        irrigate!(dFeII, z)
        # Manganese-II (solute)
        advect!(dMnII, z)
        diffuse!(dMnII, z)
        irrigate!(dMnII, z)
        # Particulate organic carbon, fast-slow-refractory (solid)
        advect!(pfoc, z)
        diffuse!(pfoc, z)
        advect!(psoc, z)
        diffuse!(psoc, z)
        advect!(proc, z)
        diffuse!(proc, z)
        # Other particulates
        advect!(pFeOH3, z)
        diffuse!(pFeOH3, z)
        advect!(pMnO2, z)
        diffuse!(pMnO2, z)
        # --- Then do the reactions! -------------------------------------------
        (
            rate_dO2,
            rate_dtCO2,
            rate_dtNO3,
            rate_dtSO4,
            rate_dtPO4,
            rate_dtNH4,
            rate_dtH2S,
            rate_dFeII,
            rate_dMnII,
            rate_pfoc,
            rate_psoc,
            rate_pFeOH3,
            rate_pMnO2,
        ) = React.rates(
            dO2.then[z],
            dtNO3.then[z],
            pMnO2.then[z],
            pFeOH3.then[z],
            dtSO4.then[z],
            dtNH4.then[z],
            dtH2S.then[z],
            dFeII.then[z],
            dMnII.then[z],
            pfoc.then[z] * kfast[z],
            psoc.then[z] * kslow[z],
            phiS_phi[z],
            RC,
            RN,
            RP,
        )
        react!(dO2, z, rate_dO2)
        react!(dtCO2, z, rate_dtCO2)
        react!(dtNO3, z, rate_dtNO3)
        react!(dtSO4, z, rate_dtSO4)
        react!(dtPO4, z, rate_dtPO4)
        react!(dtNH4, z, rate_dtNH4)
        react!(dtH2S, z, rate_dtH2S)
        react!(dFeII, z, rate_dFeII)
        react!(dMnII, z, rate_dMnII)
        react!(pfoc, z, rate_pfoc)
        react!(psoc, z, rate_psoc)
        # react!(proc, z, 0.0)  # "refractory" means it doesn't react!
        react!(pFeOH3, z, rate_pFeOH3)
        react!(pMnO2, z, rate_pMnO2)
    # ~~~ END SEDIMENT PROCESSING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Save output if we are at a savepoint
        if tsave
            dO2.save[z-1, sp+1] = dO2.now[z]
            dtCO2.save[z-1, sp+1] = dtCO2.now[z]
            dtNO3.save[z-1, sp+1] = dtNO3.now[z]
            dtSO4.save[z-1, sp+1] = dtSO4.now[z]
            dtPO4.save[z-1, sp+1] = dtPO4.now[z]
            dtNH4.save[z-1, sp+1] = dtNH4.now[z]
            dtH2S.save[z-1, sp+1] = dtH2S.now[z]
            dFeII.save[z-1, sp+1] = dFeII.now[z]
            dMnII.save[z-1, sp+1] = dMnII.now[z]
            pfoc.save[z-1, sp+1] = pfoc.now[z]
            psoc.save[z-1, sp+1] = psoc.now[z]
            proc.save[z-1, sp+1] = proc.now[z]
            pFeOH3.save[z-1, sp+1] = pFeOH3.now[z]
            pMnO2.save[z-1, sp+1] = pMnO2.now[z]
            if z == ndepths-1
                println("Radi reached savepoint $sp (step $t of $ntps)...")
                sp += 1
            end
        end
    end  # for z in 2:(ndepths-1)
    # Copy results into "previous step" arrays
    @simd for z in 2:(ndepths-1)
        dO2.then[z] = dO2.now[z]
        dtCO2.then[z] = dtCO2.now[z]
        dtNO3.then[z] = dtNO3.now[z]
        dtSO4.then[z] = dtSO4.now[z]
        dtPO4.then[z] = dtPO4.now[z]
        dtNH4.then[z] = dtNH4.now[z]
        dtH2S.then[z] = dtH2S.now[z]
        dFeII.then[z] = dFeII.now[z]
        dMnII.then[z] = dMnII.now[z]
        pfoc.then[z] = pfoc.now[z]
        psoc.then[z] = psoc.now[z]
        proc.then[z] = proc.now[z]
        pFeOH3.then[z] = pFeOH3.now[z]
        pMnO2.then[z] = pMnO2.now[z]
    end  # for z in 2:(ndepths-1)
end  # for t, main Radi model loop
# ===== End of main model loop =================================================
println("Radi done!")
# Reformat times and depths for plotting output
savetimes = [0.0 timesteps[savepoints]...]
depths_out = depths[2:end-1]
depths_out = reshape(depths_out, (length(depths_out), 1))
return (
    savetimes,
    depths_out,
    dO2.save,
    dtCO2.save,
    dtNO3.save,
    dtSO4.save,
    dtPO4.save,
    dtNH4.save,
    dtH2S.save,
    dFeII.save,
    dMnII.save,
    pfoc.save,
    psoc.save,
    proc.save,
    pFeOH3.save,
    pMnO2.save,
)
end  # function model


# "Calculate how far from equilibrium the sediment column is."
# function disequilibrium(dO2, pfoc)
#     return dO2, pfoc
# end  # function disequilibrium


sayhello() = print(raw"""

      ██▀███   ▄▄▄      ▓█████▄  ██▓
      ▓██ ▒ ██▒▒████▄    ▒██▀ ██▌▓██▒
      ▓██ ░▄█ ▒▒██  ▀█▄  ░██   █▌▒██▒
      ▒██▀▀█▄  ░██▄▄▄▄██ ░▓█▄   ▌░██░
      ░██▓ ▒██▒ ▓█   ▓██▒░▒████▓ ░██░
      ░ ▒▓ ░▒▓░ ▒▒   ▓▒█░ ▒▒▓  ▒ ░▓
       ░▒ ░ ▒░  ▒   ▒▒ ░ ░ ▒  ▒  ▒ ░
       ░░   ░   ░   ▒    ░ ░  ░  ▒ ░
        ░           ░  ░   ░     ░
                         ░

     “What do I know of man's destiny?
   I could tell you more about radishes.”
                  -- Samuel Beckett

""")

end  # module Radi
