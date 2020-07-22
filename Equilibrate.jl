module Equilibrate

const ln10 = log(10.0)

function K_NH3_CW95(TempK, Sal, SWStoTOT)
    # Ammonia dissociation constant from Clegg and Whitfield (1995)
    # via jonathansharp/CO2-System-Extd@master [2020-07-21]
    pK_NH3 = @. (
        9.244605 - 2729.33*(1.0/298.15 - 1.0/TempK + (0.04203362 - 11.24742/TempK)*Sal^0.25 +
        (-13.6416 + 1.176949*TempK^0.5 - 0.02860785TempK + 545.4834/TempK)*Sal^0.5 +
        (-0.1462507 + 0.0090226468*TempK^0.5 - 0.0001471361TempK + 10.5425/TempK)*Sal^1.5 +
        (0.004669309 - 0.0001691742*TempK^0.5 - 0.5677934/TempK)*Sal^2 +
        (-2.354039e-05 + 0.009698623/TempK)*Sal^2.5)
    )
    K_NH3 = @. 10.0 ^ -pK_NH3  # Total scale, mol/kg-H2O
    K_NH3 = @. K_NH3 * (1.0 - 0.001005Sal)  # mol/kg-SW
    K_NH3 = K_NH3 ./ SWStoTOT  # converts to SWS pH scale
end  # function K_NH3_CW95

function K_H2S_M88(TempK, Sal, SWStoTOT)
    # First hydrogen sulfide dissociation constant from Millero et al. (1988)
    # via jonathansharp/CO2-System-Extd@master [2020-07-21]
    K_H2S = @. (
        exp(225.838 - 13275.3/TempK - 34.6435*log(TempK) + 0.3449*Sal^0.5 - 0.0274*Sal)
        / SWStoTOT
    )
end  # function K_H2S_M88

function alk_borate(h::Float64, TB::Float64, KB::Float64)
    TB * KB / (KB + h)
end  # function alk_borate

function alk_silicate(h::Float64, TSi::Float64, KSi::Float64)
    TSi * KSi / (KSi + h)
end  # function alk_silicate

function alk_water(h::Float64, Kw::Float64)
    Kw / h - h
end  # function alk_water

function alk_phosphate(h::Float64, TP::Float64, KP1::Float64, KP2::Float64, KP3::Float64)
    (
        TP
        * (KP1 * KP2 * h + 2.0 * KP1 * KP2 * KP3 - h^3)
        / (
            h^3
            + KP1 * h^2
            + KP1 * KP2 * h
            + KP1 * KP2 * KP3
        )
    )
end  # function alk_phosphate

function alk_ammonia(h::Float64, TNH3::Float64, KNH3::Float64)
    TNH3 * KNH3 / (KNH3 + h)
end  # function alk_ammonia

function alk_sulfide(h::Float64, TH2S::Float64, KH2S::Float64)
    TH2S * KH2S / (KH2S + h)
end  # function alk_sulfide

function alk_sulfate(h::Float64, TSO4::Float64, KSO4::Float64)
    -TSO4 / (1.0 + KSO4 / h)
end  # function alk_sulfate

function alk_fluoride(h::Float64, TF::Float64, KF::Float64)
    -TF / (1.0 + KF / h)
end  # function alk_fluoride

function alk_carbonate(h::Float64, TC::Float64, K1::Float64, K2::Float64)
    TC*K1 * (h + 2.0K2) / (h^2 + K1*h + K1*K2)
end  # function alk_carbonate

function dalk_dpH(
    h::Float64,
    TC::Float64,
    borate_alk::Float64,
    water_alk::Float64,
    K1::Float64,
    K2::Float64,
    KB::Float64,
)
    ln10 * (TC * K1*h * (h^2 + K1*K2 + 4h*K2) / (h^2 + K1*h + K1*K2)^2 +
        borate_alk*h / (KB + h) + water_alk + 2.0h)
end  # function dalk_dpH

function dalk_dh(
    h::Float64,
    TC::Float64,
    borate_alk::Float64,
    K1::Float64,
    K2::Float64,
    KB::Float64,
    Kw::Float64,
)
    -(TC * K1 * (h^2 + K1*K2 + 4h*K2) / (h^2 + K1*h + K1*K2)^2 +
        borate_alk / (KB + h) + Kw / h^2 + 1.0)
end  # function dalk_dh

end  # module Equilibrate
