module RADI

using Base.SimdLoop
include("gsw_rho.jl")
include("Params.jl")

"Define Solute type."
struct Solute
    now::Array{Float64,1}
    then::Array{Float64,1}
    above::Float64
    dvar::Array{Float64,1}
    save::Array{Float64,2}
end # struct Solute

"Define Solid type."
struct Solid
    now::Array{Float64,1}
    then::Array{Float64,1}
    dvar::Array{Float64,1}
    save::Array{Float64,2}
end # struct Solute

SoluteOrSolid = SolidOrSolute = Union{Solid,Solute}
FloatOrArray = ArrayOrFloat = Union{Float64,Array{Float64,1}}

"Constructor for a Solute."
function Solute(var_start::Array{Float64,1}, above::Float64,
        dvar::Array{Float64,1}, var_save::Array{Float64,2})
    return Solute(copy(var_start), copy(var_start), above, dvar, var_save)
end # function Solute

"Constructor for a Solid."
function Solid(var_start::Array{Float64,1}, dvar::Array{Float64,1},
        var_save::Array{Float64,2})
    return Solid(copy(var_start), copy(var_start), dvar, var_save)
end # function Solid

"Prepare vectors of model timesteps and savepoints. All time units are in days."
function preptime(stoptime::Float64, interval::Float64, saveperXsteps::Int)
    timesteps::Array{Float64,1} = collect(0.0:interval:stoptime)
    ntps::Int = length(timesteps)
    savepoints::Array{Int64,1} = collect(1:saveperXsteps:ntps)
    # Save final timepoint if it's not already in the list
    if !(ntps in savepoints)
        append!(savepoints, ntps)
    end
    nsps::Int = length(savepoints)
    return timesteps, savepoints, ntps, nsps
end # function preptime

"Prepare the model's depth vector. All depth units are in metres."
function prepdepth(depth_res::Float64, depth_max::Float64)
    depths = collect(-depth_res:depth_res:(depth_max+depth_res)) # in m
    depths[1] = NaN
    depths[end] = NaN
    ndepths = length(depths)
    depth_res2 = depth_res^2
    return depths, ndepths, depth_res2
end # function prepdepth

"Assemble all RADI model parameters."
function assemble(
    depths::Array{Float64},
    z_res::Float64,
    dtPO4_w::Float64,
    dO2_w::Float64,
    T::Float64,
    S::Float64,
    P::Float64,
    dbl::Float64,
    Fpom::Float64,
    rho_pom::Float64,
    phi0::Float64,
    phiInf::Float64,
    beta::Float64,
    lambda_b::Float64,
    lambda_i::Float64,
)

    # "Redfield" ratios and OM stoichiometry
    rho_sw = gsw_rho(S, T, P) # seawater density [kg/m^3]
    RC, RN, RP = Params.redfield(dtPO4_w, rho_sw)
    Mpom = Params.rmm_pom(RC, RN, RP)
    Fpom_mol = Fpom/Mpom
    Fpoc = Fpom_mol*RC

    # Depth-dependent porosity
    phi = Params.phi(phi0, phiInf, beta, depths)
    phiS = Params.phiS(phi) # solid volume fraction
    phiS_phi = phiS./phi
    tort2 = Params.tort2(phi) # tortuosity squared from Boudreau (1996, GCA)
    delta_phi = Params.delta_phi(phi0, phiInf, beta, depths)
    delta_phiS = Params.delta_phiS(delta_phi)
    delta_tort2i = Params.delta_tort2i(delta_phi, phi, tort2)
    delta_tort2i_tort2 = delta_tort2i.*tort2

    # Bioturbation (for solids)
    D_bio_0 = Params.D_bio_0(Fpoc)
    # ^[m2/a] surf bioturb coeff, Archer et al. (2002)
    D_bio = Params.D_bio(depths, D_bio_0, lambda_b, dO2_w)
    # ^[m2/a] bioturb coeff, Archer et al (2002)
    delta_D_bio = Params.delta_D_bio(depths, D_bio, lambda_b)

    # Organic matter degradation parameters
    krefractory = Params.krefractory(depths, D_bio_0)
    # ^[/a] from Archer et al (2002)

    # Solid fluxes and solid initial conditions
    x0 = Params.x0(Fpom, rho_pom, phiS[2])
    # ^[m/a] bulk burial velocity at sediment-water interface
    xinf = Params.xinf(x0, phiS[2], phiS[end-1])
    # ^[m/a] bulk burial velocity at the infinite depth
    u = Params.u(xinf, phi) # [m/a] porewater burial velocity
    w = Params.w(xinf, phiS) # [m/a] solid burial velocity

    # Biodiffusion depth-attenuation: see Boudreau (1996); Fiadeiro and Veronis
    # (1977)
    Peh = Params.Peh(w, z_res, D_bio)
    # ^one half the cell Peclet number (Eq. 97 in Boudreau 1996)
    # when Peh<<1, biodiffusion dominates, when Peh>>1, advection dominates
    sigma = Params.sigma(Peh)
    sigma1m = 1.0 .- sigma
    sigma1p = 1.0 .+ sigma

    # vvv NOT YET IN THE PARAMETERS PART OF THE DOCUMENTATION vvvvvvvvvvvvvvvvvv
    # Temperature-dependent "free solution" diffusion coefficients
    D_dO2 = Params.D_dO2(T)
    # ^[m2/a] oxygen diffusion coefficient from Li and Gregory
    D_dO2_tort2 = D_dO2./tort2

    # Irrigation (for solutes)
    alpha_0 = Params.alpha_0(Fpoc, dO2_w)
    # ^[/a] from Archer et al (2002)
    alpha = Params.alpha(alpha_0, depths, lambda_i)
    # ^[/a] Archer et al (2002)

    APPW = Params.APPW(w, delta_D_bio, delta_phiS, D_bio, phiS)
    TR = Params.TR(z_res, tort2[2], dbl)
    Fpoc_phiS_0 = Fpoc/phiS[2]
    zr_Db_0 = 2.0z_res/D_bio[2]
    # ^^^ NOT YET IN THE PARAMETERS PART OF THE DOCUMENTATION ^^^^^^^^^^^^^^^^^^

    return Fpoc, phi, phiS, phiS_phi, tort2, delta_phi, delta_phiS,
        delta_tort2i_tort2, D_bio, krefractory, u, w, sigma, sigma1m, sigma1p,
        D_dO2_tort2, alpha, APPW, TR, Fpoc_phiS_0, zr_Db_0

end # function assemble

"Run the RADI model"
function model(
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
    lambda_i::Float64,
    T::Float64,
    S::Float64,
    P::Float64,
    dO2_w::Float64,
    dtPO4_w::Float64,
    Fpom::Float64,
    rho_pom::Float64,
    dO2_i::FloatOrArray,
    poc_i::FloatOrArray,
)

# Set up model time and depth grids
timesteps, savepoints, ntps, nsps = preptime(stoptime, interval, saveperXsteps)
depths, ndepths, z_res2 = prepdepth(z_res, z_max)
sp = 1 # initialise savepoints

# Assemble model parameters
Fpoc, phi, phiS, phiS_phi, tort2, delta_phi, delta_phiS, delta_tort2i_tort2,
        D_bio, krefractory, u, w, sigma, sigma1m, sigma1p, D_dO2_tort2, alpha,
        APPW, TR, Fpoc_phiS_0, zr_Db_0 =
    assemble(depths, z_res, dtPO4_w, dO2_w, T, S, P, dbl, Fpom, rho_pom,
        phi0, phiInf, beta, lambda_b, lambda_i)

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
    return then[3] + (Fpoc_phiS_0 - w[2]*then[2])*zr_Db_0
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
        delta_phi_z::Float64, phi_z::Float64, delta_tort2i_tort2_z::Float64,
        D_var::Float64)
    return -(u_z - delta_phi_z*D_var/phi_z -
        D_var*delta_tort2i_tort2_z)*(then_z1p - then_z1m)/(2.0z_res)
end # function advect

"Advect a Solute."
function advect!(z::Int, var::Solute)
    var.now[z] += interval*advectsolute(var.then[z+1], var.then[z-1], u[z],
        delta_phi[z], phi[z], delta_tort2i_tort2[z], var.dvar[z])
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
    var.now[z] += interval*-APPW[z]*(sigma1m[z]*var.then[z+1] +
        2.0sigma[z]*var.then[z] - sigma1p[z]*var.then[z-1])/(2.0z_res)
# No idea why the approach below is so much slower, only for this function?!
    # var.now[z] += interval*advectsolid(var.then[z], var.then[z+1],
    #     var.then[z-1], APPW[z], sigma[z], sigma1p[z], sigma1m[z])
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
function irrigate(then_z::Float64, above::Float64, alpha_z::Float64)
    return alpha_z*(above - then_z)
end # function irrigate

"Irrigate a Solute throughout the sediment."
function irrigate!(z::Int, var::Solute)
    var.now[z] += interval*irrigate(var.then[z], var.above, alpha[z])
end # function irrigate!

# ===== Run RADI run! ==========================================================
# Create variables to model
dO2 = makeSolute(dO2_i, dO2_w, D_dO2_tort2)
poc = makeSolid(poc_i, D_bio)
# Main RADI model loop
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
        # --- Then do the reactions! -------------------------------------------
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

"Calculate how far from equilibrium the sediment column is."
function disequilibrium(dO2, poc)
    return dO2, poc
end # function disequilibrium

end # module RADI
