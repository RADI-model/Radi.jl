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

    timesteps, savepoints, ntps, nsps = timeprep(stoptime, interval, saveperXsteps)
    sp = 1

    # Set water column conditions
    T = 1.4 # temperature / degC
    S = 34.69 # practical salinity
    rho_sw = gsw_rho(S, T, 1) # seawater density [kg/m^3]
    oxy_w = 159.7e-6*rho_sw # dissolved oxygen [mol/m3]
    po4_w = 2.39e-6*rho_sw # phosphate from GLODAP at station location, bottom
    # waters [mol/m3]

    # "Redfield" ratios
    RC = @. 1/(6.9e-3po4_w/(1e-6rho_sw) + 6e-3)
    # ^P:C computed as a function of SRP from Galbraith and Martiny PNAS 2015
    RN = 11 # value at 60 degS from Martiny et al. Nat G 2013
    RP = 1 # Redfield ratio for P in the deep sea
    M_OM = @. 30.03 + (RN/RC)*17.03 + (RP/RC)*97.994
    # ^ Organic Matter molar mass [g of OM per mol of OC]

    # Create arrays for modelled variables
    z_res = 0.001
    z_res2 = z_res^2
    depths = 0:z_res:0.2 # in m
    ndepths = length(depths)
    oxy0 = fill(oxy_w, (ndepths,))
    oxy = copy(oxy0)
    poc0 = fill(0.0001, (ndepths,))
    poc = copy(poc0)

    # Depth-dependent porosity
    phi = @. (0.85 - 0.74)*exp(-33*depths) + 0.74 # porosity profile (porewater
    # bulk fraction) fitted from station7 mooring3 of cruise NBP98-2 by Sayles
    # et al. DSR 2001
    phiS = 1 .- phi # solid volume fraction
    phiS_phi = phiS ./ phi
    tort = @. sqrt(1 - 2log(phi)) # tortuosity from Boudreau (1996, GCA)
    tort2 = tort.^2 # tortuosity squared

    # Solid fluxes and solid initial conditions
    Foc = 1.0 # flux of total organic carbon to the bottom [mol/m2/a]
    Ftot = Foc.*M_OM # total sediment flux [g/m2/a]
    v0 = Ftot/(2.65e6*phiS[1]) # [m/a] bulk burial velocity at sediment-water interface
    vinf = v0*phiS[1]/phiS[end] # [m/a] bulk burial velocity at the infinite depth
    u = @. vinf*phi[end]/phi # [m/a] porewater burial velocity
    w = @. vinf*phiS[end]/phiS # [m/a] solid burial velocity

    # Diffusive boundary layer
    dbl = 1e-3 # thickness at location taken from Sulpis et al. 2018 PNAS [m]

    # temperature dependent "free solution" diffusion coefficients
    D_O2 = @. 0.034862 + 0.001409T # [m2/a] oxygen diffusion coefficient from Li and Gregiry

    # bioturbation (for solids)
    D_bio_0 = @. 0.0232e-4*(Foc*1e2)^0.85 # [m2/a] surf bioturb coeff, Archer et al (2002)
    D_bio = @. D_bio_0*exp(-(depths/0.08)^2)*((oxy_w/1e-3)/((oxy_w/1e-3)+20))
    # ^[m2/a] bioturb coeff, Archer et al (2002)

    # irrigation (for solutes)
    alpha_0 = @. 11*(atan((5*Foc*1e2 - 400)/400)/3.1416 + 0.5) - 0.9 +
        20*((oxy_w/1e-3)/((oxy_w/1e-3) + 10))*exp(-(oxy_w/1e-3)/10)*Foc*1e2/(Foc*1e2+30)
    # ^[/a] from Archer et al (2002)
    alpha = @. alpha_0*exp(-(depths/0.05)^2)
    # ^[/a] Archer et al (2002) the depth of 5 cm was changed

    # Depth-dependent porosity and diffusion coefficient loss
    delta_phi = append!([0.0], diff(phi)) # depth-dependent porosity loss
    delta_phiS = append!([0.0], diff(phiS)) # depth-dependent solid fraction gain
    delta_tort2 = append!([0.0], diff(tort.^2)) # depth-dependent tortuosity gain
    delta_D_bio = append!([0.0], diff(D_bio)) # [m/a]

    # Biodiffusion depth-attenuation: see Boudreau (1996); Fiadeiro and Veronis (1977)
    Peh = @. w*z_res/(2D_bio) # one half the cell Peclet number (Eq. 97 in Boudreau 1996)
    # when Peh<<1, biodiffusion dominates, when Peh>>1, advection dominates
    sigma = @. 1/tanh(Peh) - 1/(Peh) # Eq. 96 in Boudreau 1996

    # Organic matter degradation parameters
    krefractory = @. 80.25D_bio_0*exp(-depths) # [/a] from Archer et al (2002)

    # Short-cut transport variables
    APPW = @. w - delta_D_bio - delta_phiS*D_bio/phiS
    DFF = @. (tort2*delta_phi/phi - delta_tort2)/tort2^2
    DBF = @. D_bio*(phiS - delta_phi)
    TR = @. 2z_res*tort^2/dbl

    # Create arrays for saved variables
    oxy_save = fill(NaN, (ndepths, nsps))
    poc_save = fill(NaN, (ndepths, nsps))

    function react!(z::Int, var::Array{Float64,1}, rate::Float64)
        var[z] += interval*rate
    end # function react!

    function diffuse!(z::Int, var0::Array{Float64,1}, var::Array{Float64,1})
        var[z] += interval*(var0[z+1] - 2var0[z] + var0[z-1])*
            2D_O2/(tort2[z]*z_res2)
    end # function diffuse!

    # ===== Run RADI run: main model loop ======================================
    for t in 1:ntps
        tsave = t in savepoints # do we save this time?

        for z in 1:ndepths
            # Redox reaction rates
            Rg_z = poc[z]*krefractory[z]

            # Reactions
            react!(z, oxy, -Rg_z*phiS_phi[z])
            react!(z, poc, -Rg_z)

            # Diffusion
            if z > 1 && z < ndepths
                diffuse!(z, oxy0, oxy)
            end # if

            # Save output if we are at a savepoint
            if tsave
                if t==1
                    oxy_save[z, sp] = oxy0[z]
                    poc_save[z, sp] = poc0[z]
                else
                    oxy_save[z, sp] = oxy[z]
                    poc_save[z, sp] = poc[z]
                end
                if z==ndepths
                    println("RADI: reached savepoint $sp (step $t of $ntps)...")
                    sp += 1
                end # if
            end # if

            # Save result to "previous step" arrays
            oxy0[z] = oxy[z]
            poc0[z] = poc[z]

        end # for z
    end # for t

    return depths, oxy_save, poc_save

end # function model

say_RADI() = println("RADI done!")

end # module RADI
