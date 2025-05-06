module React

# Monod scheme constants in mol/m^3
# For oxygen and nitrate from Soetaert et al. 1996 (GCA)
const KM_dO2 = 0.003
const KMi_dO2 = 0.01
const KM_dtNO3 = 0.03
const KMi_dtNO3 = 0.005
# For others from van Cappellen and Wang (1996)
const KM_pMnO2 = 42.4
const KMi_pMnO2 = 42.4
const KM_pFeOH3 = 265.0
const KMi_pFeOH3 = 265.0
const KM_dtSO4 = 1.6
const KMi_dtSO4 = 1.6

# Monod scheme inhibition factors
inhibition_dO2(dO2::Float64) = KMi_dO2 / (KMi_dO2 + dO2)
inhibition_dtNO3(dtNO3::Float64) = KMi_dtNO3 / (KMi_dtNO3 + dtNO3)
inhibition_pMnO2(pMnO2::Float64) = KMi_pMnO2 / (KMi_pMnO2 + pMnO2)
inhibition_pFeOH3(pFeOH3::Float64) = KMi_pFeOH3 / (KMi_pFeOH3 + pFeOH3)
inhibition_dtSO4(dtSO4::Float64) = KMi_dtSO4 / (KMi_dtSO4 + dtSO4)


"""Organic matter degradation pathway factors.
From the code of Couture et al. (EST 2010), following Boudreau (1996).
"""
function degradationfactors(
    dO2::Float64,
    dtNO3::Float64,
    pMnO2::Float64,
    pFeOH3::Float64,
    dtSO4::Float64,
)
    # Evaluate inhibition factors
    Mi_dO2 = inhibition_dO2(dO2)
    Mi_dtNO3 = inhibition_dtNO3(dtNO3)
    Mi_pMnO2 = inhibition_pMnO2(pMnO2)
    Mi_pFeOH3 = inhibition_pFeOH3(pFeOH3)
    Mi_dtSO4 = inhibition_dtSO4(dtSO4)
    # Calculate degradation pathway factors
    fdO2 = dO2 / (KM_dO2 + dO2)
    fdtNO3 = Mi_dO2 *
        dtNO3 / (KM_dtNO3 + dtNO3)
    fpMnO2 = Mi_dO2 * Mi_dtNO3 *
        pMnO2 / (KM_pMnO2 + pMnO2)
    fpFeOH3 = Mi_dO2 * Mi_dtNO3 * Mi_pMnO2 *
        pFeOH3 / (KM_pFeOH3 + pFeOH3)
    fdtSO4 = Mi_dO2 * Mi_dtNO3 * Mi_pMnO2 * Mi_pFeOH3 *
        dtSO4 / (KM_dtSO4 + dtSO4)
    fdCH4 = Mi_dO2 * Mi_dtNO3 * Mi_pMnO2 * Mi_pFeOH3 * Mi_dtSO4
    fox = fdO2 + fdtNO3 + fpMnO2 + fpFeOH3 + fdtSO4 + fdCH4

    # Normalize factors
    fdO2N = fdO2 /= fox
    fdtNO3N = fdtNO3  /= fox
    fpMnO2N = fpMnO2 /= fox
    fpFeOH3N = fpFeOH3 /= fox
    fdtSO4N = fdtSO4 /= fox
    fdCH4N = fdCH4 /= fox


    foxN = fdO2N + fdtNO3N + fpMnO2N + fpFeOH3N + fdtSO4N + fdCH4N



    return fdO2N, fdtNO3N, fpMnO2N, fpFeOH3N, fdtSO4N, fdCH4N, foxN
end # function degradationfactors


function degrade(
    dO2::Float64,
    dtNO3::Float64,
    pMnO2::Float64,
    pFeOH3::Float64,
    dtSO4::Float64,
    pfoc_kfast::Float64,
    psoc_kslow::Float64,
)
    # Evaluate degradation factors
    fdO2, fdtNO3, fpMnO2, fpFeOH3, fdtSO4, fdCH4, fox =
        degradationfactors(dO2, dtNO3, pMnO2, pFeOH3, dtSO4)
    # Calculate degradation reaction rates
    Rfast_dO2 = fdO2 * pfoc_kfast
    Rslow_dO2 = fdO2 * psoc_kslow
    Rfast_dtNO3 = fdtNO3 * pfoc_kfast
    Rslow_dtNO3 = fdtNO3 * psoc_kslow
    Rfast_pMnO2 = fpMnO2 * pfoc_kfast
    Rslow_pMnO2 = fpMnO2 * psoc_kslow
    Rfast_pFeOH3 = fpFeOH3 * pfoc_kfast
    Rslow_pFeOH3 = fpFeOH3 * psoc_kslow
    Rfast_dtSO4 = fdtSO4 * pfoc_kfast
    Rslow_dtSO4 = fdtSO4 * psoc_kslow
    Rfast_dCH4 = fdCH4 * pfoc_kfast
    Rslow_dCH4 = fdCH4 * psoc_kslow
    Rfast_total = fox * pfoc_kfast
    Rslow_total = fox * psoc_kslow
    return (
        Rfast_dO2,
        Rslow_dO2,
        Rfast_dtNO3,
        Rslow_dtNO3,
        Rfast_pMnO2,
        Rslow_pMnO2,
        Rfast_pFeOH3,
        Rslow_pFeOH3,
        Rfast_dtSO4,
        Rslow_dtSO4,
        Rfast_dCH4,
        Rslow_dCH4,
        Rfast_total,
        Rslow_total,
    )
end # function degrade


# Redox reaction first order rate constants for deep sea from Boudreau (1996)
# All in mol/m^3/a
const kMnII_redox = 1e6
const kFeII_redox = 1e6
const kNH3_redox = 1e4
const kH2S_redox = 3e5
const kCH4_O2redox = 1e7
const kCH4_SO4redox = 1e7

"Redox reaction rates."
function redox(
    dO2::Float64,
    dtNH4::Float64,
    dtH2S::Float64,
    dFeII::Float64,
    dMnII::Float64,
    dCH4::Float64,
    dtSO4::Float64,
)
    R_NH3_redox = kNH3_redox * dtNH4 * dO2
    R_H2S_redox = kH2S_redox * dtH2S * dO2
    R_MnII_redox = kMnII_redox * dMnII * dO2
    R_FeII_redox = kFeII_redox * dFeII * dO2
    R_CH4_O2redox = kCH4_O2redox * dCH4 * dO2
    R_CH4_SO4redox = kCH4_SO4redox * dCH4 * dtSO4
    return R_MnII_redox, R_FeII_redox, R_NH3_redox, R_H2S_redox, R_CH4_O2redox, R_CH4_SO4redox
end # function redox


"CaCO3 mineral dissolution and precipitation rates."
function dissolve_precipitate_CaCO3(
    pcalcite::Float64,
    paragonite::Float64,
    dCa::Float64,
    dCO3::Float64,
    KCa::Float64,
    KAr::Float64,
)
    # Saturation states
    OmegaCa = dCa * dCO3 / KCa
    OmegaAr = dCa * dCO3 / KAr

    # Calcite dissolution rate from Naviaux et al. (2019) Marine Chemistry
    if 0.8275 < OmegaCa <= 1.0
        Rdiss_calcite = pcalcite * 0.00632 * (1.0 - OmegaCa) ^ 0.11
    elseif OmegaCa <= 0.8275  # this is the Omega value for which both laws are equal
        Rdiss_calcite = pcalcite * 20.0 * (1.0 - OmegaCa) ^ 4.7
    else
        Rdiss_calcite = 0.0
    end  # calcite dissolution rate
    # Aragonite dissolution rate from Dong et al. (2019) EPSL
    if 0.835 < OmegaAr <= 1.0
        Rdiss_aragonite = paragonite * 0.0038000000000000004 * (1.0 - OmegaAr) ^ 0.13
    elseif OmegaAr <= 0.835
        Rdiss_aragonite = paragonite * 0.042 * (1.0 - OmegaAr) ^ 1.46
    else
        Rdiss_aragonite = 0.0
    end   # aragonite dissolution rate
    # Calcite precipitation rate from Zuddas and Mucci, GCA (1998),
    # normalised to the same surface area as dissolution (4m^2/g)
    if OmegaCa > 1.0
        Rprec_calcite = 0.4075 * (OmegaCa - 1.0) ^ 1.76
    else
        Rprec_calcite = 0.0
    end
    # Aragonite does not currently precipitate
    Rprec_aragonite = 0.0
    return Rdiss_calcite, Rdiss_aragonite, Rprec_calcite, Rprec_aragonite
end  # function dissolve_precipitate_CaCO3

"All reaction rates."
function getreactions(
    dO2::Float64,
    dtNO3::Float64,
    pMnO2::Float64,
    pFeOH3::Float64,
    dtSO4::Float64,
    dtNH4::Float64,
    dtH2S::Float64,
    dFeII::Float64,
    dMnII::Float64,
    dCH4::Float64,
    pfoc_kfast::Float64,
    psoc_kslow::Float64,
    pcalcite::Float64,
    paragonite::Float64,
    dCa::Float64,
    dCO3::Float64,
    KCa::Float64,
    KAr::Float64,
)
    (
        Rfast_dO2,
        Rslow_dO2,
        Rfast_dtNO3,
        Rslow_dtNO3,
        Rfast_pMnO2,
        Rslow_pMnO2,
        Rfast_pFeOH3,
        Rslow_pFeOH3,
        Rfast_dtSO4,
        Rslow_dtSO4,
        Rfast_dCH4,
        Rslow_dCH4,
        Rfast_total,
        Rslow_total,
    ) = degrade(dO2, dtNO3, pMnO2, pFeOH3, dtSO4, pfoc_kfast, psoc_kslow)
    # Redox reactions
    R_dMnII, R_dFeII, R_dNH3, R_dH2S, R_CH4_O2redox,  R_CH4_SO4redox = redox(dO2, dtNH4, dtH2S, dFeII, dMnII, dCH4, dtSO4)
    Rdiss_calcite, Rdiss_aragonite, Rprec_calcite, Rprec_aragonite = dissolve_precipitate_CaCO3(pcalcite, paragonite, dCa, dCO3, KCa, KAr)
    return (
        Rfast_dO2,
        Rslow_dO2,
        Rfast_dtNO3,
        Rslow_dtNO3,
        Rfast_pMnO2,
        Rslow_pMnO2,
        Rfast_pFeOH3,
        Rslow_pFeOH3,
        Rfast_dtSO4,
        Rslow_dtSO4,
        Rfast_dCH4,
        Rslow_dCH4,
        Rfast_total,
        Rslow_total,
        R_dMnII,
        R_dFeII,
        R_dNH3,
        R_dH2S,
        R_CH4_O2redox,
        R_CH4_SO4redox,
        Rdiss_calcite,
        Rdiss_aragonite,
        Rprec_calcite,
        Rprec_aragonite,
    )
end # function getreactions


"Convert reactions to individual component rates."
function reactions2rates(
    Rfast_dO2::Float64,
    Rslow_dO2::Float64,
    Rfast_dtNO3::Float64,
    Rslow_dtNO3::Float64,
    Rfast_pMnO2::Float64,
    Rslow_pMnO2::Float64,
    Rfast_pFeOH3::Float64,
    Rslow_pFeOH3::Float64,
    Rfast_dtSO4::Float64,
    Rslow_dtSO4::Float64,
    Rfast_dCH4::Float64,
    Rslow_dCH4::Float64,
    Rfast_total::Float64,
    Rslow_total::Float64,
    R_dMnII::Float64,
    R_dFeII::Float64,
    R_dNH3::Float64,
    R_dH2S::Float64,
    R_CH4_O2redox::Float64,
    R_CH4_SO4redox::Float64,
    Rdiss_calcite::Float64,
    Rdiss_aragonite::Float64,
    Rprec_calcite::Float64,
    Rprec_aragonite::Float64,
    phiS_phi_z::Float64,
    RC::Float64,
    RN::Float64,
    RP::Float64,
)
    # Add things up for convenience
    Rdeg_dO2 = Rfast_dO2 + Rslow_dO2
    Rdeg_dtNO3 = Rfast_dtNO3 + Rslow_dtNO3
    Rdeg_dtSO4 = Rfast_dtSO4 + Rslow_dtSO4
    Rdeg_pFeOH3 = Rfast_pFeOH3 + Rslow_pFeOH3
    Rdeg_pMnO2 = Rfast_pMnO2 + Rslow_pMnO2
    Rdeg_dCH4 = Rfast_dCH4 + Rslow_dCH4
    Rdeg_total = Rfast_total + Rslow_total
    # Net CaCO3 dissolution
    Rdiss_CaCO3 = Rdiss_calcite + Rdiss_aragonite - Rprec_calcite - Rprec_aragonite
    # Total changes in porewater/sediment components from reaction rates
    p2d = phiS_phi_z  # convert particulate to dissolved
    d2p = 1.0 / phiS_phi_z  # convert dissolved to particulate
    rate_dO2 = p2d * -Rdeg_dO2 - (R_dFeII/4.0 + R_dMnII/2.0 + R_dH2S*2.0 + R_dNH3*2.0) - 2.0*R_CH4_O2redox
    rate_dtCO2 = p2d * (RC * (Rdeg_total - Rdeg_dCH4/2.0) + Rdiss_CaCO3) + R_CH4_O2redox + R_CH4_SO4redox
    rate_dtNO3 = p2d * RC * -Rdeg_dtNO3*0.8 + R_dNH3
    rate_dtSO4 = p2d * RC * -Rdeg_dtSO4/2.0 + R_dH2S - R_CH4_SO4redox
    rate_dCH4 = p2d * RC * Rdeg_dCH4/2.0 - R_CH4_O2redox - R_CH4_SO4redox
    rate_dtPO4 = p2d * RP * Rdeg_total
    rate_dtNH4 = p2d * RN * Rdeg_total - R_dNH3
    rate_dtH2S = p2d * RC * Rdeg_dtSO4/2.0 - R_dH2S + R_CH4_SO4redox
    rate_dFeII = p2d * RC * Rdeg_pFeOH3*4.0 - R_dFeII
    rate_dMnII = p2d * RC * Rdeg_pMnO2*2.0 - R_dMnII
    rate_pfoc = -Rfast_total
    rate_psoc = -Rslow_total
    rate_pFeOH3 = RC * -Rdeg_pFeOH3*4.0 + d2p * R_dFeII
    rate_pMnO2 = RC * -Rdeg_pMnO2*2.0 + d2p * R_dMnII
    # Total alkalinity and calcium
    rate_dCa = p2d * Rdiss_CaCO3
    rate_dalk = p2d * (
            (RN - RP) * (Rdeg_dO2 + Rdeg_dCH4) +
            (RN - RP + 0.8RC) * Rdeg_dtNO3 +
            (RN - RP + 4.0RC) * Rdeg_pMnO2 +
            (RN - RP + 8.0RC) * Rdeg_pFeOH3 +
            (RN - RP + 1.0RC) * Rdeg_dtSO4 +
            2.0 * Rdiss_CaCO3) + 
            2.0 * R_CH4_SO4redox -
            2.0 * (R_dMnII + R_dFeII + R_dNH3 + R_dH2S)
    # CaCO3 minerals
    rate_pcalcite = - Rdiss_calcite + Rprec_calcite
    rate_paragonite = - Rdiss_aragonite
    return (
        rate_dO2,
        rate_dtCO2,
        rate_dtNO3,
        rate_dtSO4,
        rate_dtPO4,
        rate_dtNH4,
        rate_dtH2S,
        rate_dFeII,
        rate_dMnII,
        rate_dCH4,
        rate_dalk,
        rate_dCa,
        rate_pfoc,
        rate_psoc,
        rate_pFeOH3,
        rate_pMnO2,
        rate_pcalcite,
        rate_paragonite,

        Rdeg_dO2, 
        Rdeg_dtNO3, 
        Rdeg_dtSO4, 
        Rdeg_pFeOH3, 
        Rdeg_pMnO2, 
        Rdeg_dCH4,
        Rdeg_total
    )
end # function reactions2rates


"Rates of change of each component."
function rates(
    dO2,
    dtNO3,
    pMnO2,
    pFeOH3,
    dtSO4,
    dtNH4,
    dtH2S,
    dFeII,
    dMnII,
    dCH4,
    pfoc_kfast,
    psoc_kslow,
    pcalcite,
    paragonite,
    dCa,
    dCO3,
    KCa,
    KAr,
    phiS_phi_z,
    RC,
    RN,
    RP,
)
    # Get total reaction rates from concentrations
    (
        Rfast_dO2,
        Rslow_dO2,
        Rfast_dtNO3,
        Rslow_dtNO3,
        Rfast_pMnO2,
        Rslow_pMnO2,
        Rfast_pFeOH3,
        Rslow_pFeOH3,
        Rfast_dtSO4,
        Rslow_dtSO4,
        Rfast_dCH4,
        Rslow_dCH4,
        Rfast_total,
        Rslow_total,
        R_dMnII,
        R_dFeII,
        R_dNH3,
        R_dH2S,
        R_dCH4_O2redox,
        R_dCH4_SO4redox,
        Rdiss_calcite,
        Rdiss_aragonite,
        Rprec_calcite,
        Rprec_aragonite,
    ) = getreactions(
        dO2,
        dtNO3,
        pMnO2,
        pFeOH3,
        dtSO4,
        dtNH4,
        dtH2S,
        dFeII,
        dMnII,
        dCH4,
        pfoc_kfast,
        psoc_kslow,
        pcalcite,
        paragonite,
        dCa,
        dCO3,
        KCa,
        KAr,
    )
    # Convert total reaction rates to individual component rates
    return reactions2rates(
        Rfast_dO2,
        Rslow_dO2,
        Rfast_dtNO3,
        Rslow_dtNO3,
        Rfast_pMnO2,
        Rslow_pMnO2,
        Rfast_pFeOH3,
        Rslow_pFeOH3,
        Rfast_dtSO4,
        Rslow_dtSO4,
        Rfast_dCH4,
        Rslow_dCH4,
        Rfast_total,
        Rslow_total,
        R_dMnII,
        R_dFeII,
        R_dNH3,
        R_dH2S,
        R_dCH4_O2redox,
        R_dCH4_SO4redox,
        Rdiss_calcite,
        Rdiss_aragonite,
        Rprec_calcite,
        Rprec_aragonite,
        phiS_phi_z,
        RC,
        RN,
        RP,
    )
end # function rates


end # module React
