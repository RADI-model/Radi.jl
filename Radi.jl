module Radi

using MAT
include("gsw_rho.jl")
include("Model.jl")

# Import site-specific settings
include("IC_W29.jl")

"Convenient wrapper function for running RADI."
function go(initial::Dict)
    Model.sayhello()
    @time (
        savetimes,
        depths,
        dO2,
        dtCO2,
        dtNO3,
        dtSO4,
        dtPO4,
        dtNH4,
        dtH2S,
        dFeII,
        dMnII,
        dalk,
        dCa,
        pfoc,
        psoc,
        proc,
        pFeOH3,
        pMnO2,
    ) = Model.timeloop(
        stoptime,
        interval,
        saveperXsteps,
        z_max,
        z_res,
        dbl,
        phiInf,
        phi0,
        beta,
        lambda_b,
        lambda_f,
        lambda_s,
        lambda_i,
        T,
        S,
        P,
        dO2_w,
        dtCO2_w,
        dtNO3_w,
        dtSO4_w,
        dtPO4_w,
        dtNH4_w,
        dtH2S_w,
        dFeII_w,
        dMnII_w,
        dalk_w,
        dCa_w,
        Fpom,
        Fpom_r,
        Fpom_s,
        Fpom_f,
        FMnO2,
        FFeOH3,
        rho_p,
        initial[:dO2],
        initial[:dtCO2],
        initial[:dtNO3],
        initial[:dtSO4],
        initial[:dtPO4],
        initial[:dtNH4],
        initial[:dtH2S],
        initial[:dFeII],
        initial[:dMnII],
        initial[:dalk],
        initial[:dCa],
        initial[:pfoc],
        initial[:psoc],
        initial[:proc],
        initial[:pFeOH3],
        initial[:pMnO2],
    )
    return Dict(
        :savetimes => savetimes,
        :depths => depths,
        :dO2 => dO2,
        :dtCO2 => dtCO2,
        :dtNO3 => dtNO3,
        :dtSO4 => dtSO4,
        :dtPO4 => dtPO4,
        :dtNH4 => dtNH4,
        :dtH2S => dtH2S,
        :dFeII => dFeII,
        :dMnII => dMnII,
        :dalk => dalk,
        :dCa => dCa,
        :pfoc => pfoc,
        :psoc => psoc,
        :proc => proc,
        :pFeOH3 => pFeOH3,
        :pMnO2 => pMnO2,
    )
end  # function go


"Run Radi again starting from the end of a previous result."
function again(results::Dict)
    initials = (
        :dO2,
        :dtCO2,
        :dtNO3,
        :dtSO4,
        :dtPO4,
        :dtNH4,
        :dtH2S,
        :dFeII,
        :dMnII,
        :dalk,
        :dCa,
        :pfoc,
        :psoc,
        :proc,
        :pFeOH3,
        :pMnO2,
    )
    initial = Dict(i => results[i][:, end] for i in initials)
    go(initial)
end  # function again


"Save results for plotting elsewhere with an addition to the filename."
function save(results::Dict, suffix::String)
    matwrite("results/" * modelrun * suffix * ".mat", Dict(
        string(k) => v for (k, v) in results
    ))
end  # function save


"Save results for plotting elsewhere."
function save(results::Dict)
    save(results, "")
end  # function save

    
# Run the model for the first time
initial = Dict(
    :dO2 => dO2_i,
    :dtCO2 => dtCO2_i,
    :dtNO3 => dtNO3_i,
    :dtSO4 => dtSO4_i,
    :dtPO4 => dtPO4_i,
    :dtNH4 => dtNH4_i,
    :dtH2S => dtH2S_i,
    :dFeII => dFeII_i,
    :dMnII => dMnII_i,
    :dalk => dalk_i,
    :dCa => dCa_i,
    :pfoc => pfoc_i,
    :psoc => psoc_i,
    :proc => proc_i,
    :pFeOH3 => pFeOH3_i,
    :pMnO2 => pMnO2_i,
)
results = go(initial)
save(results)


end  # module Radi
