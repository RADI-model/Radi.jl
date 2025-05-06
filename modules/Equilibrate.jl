module Equilibrate

const ln10 = log(10.0)

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
