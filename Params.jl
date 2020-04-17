module Params

"Evaluate the 'Redfield' ratios for particulate organic matter."
function redfield(dtPO4_w::Float64, rho_sw::Float64)
    RC = @. 1.0/(6.9e-3dtPO4_w/(1e-6rho_sw) + 6e-3)
    # ^P:C computed as a function of SRP from Galbraith and Martiny PNAS 2015
    RN = 11.0 # value at 60 degS from Martiny et al. Nat G 2013
    RP = 1.0 # Redfield ratio for P in the deep sea
    return RC, RN, RP
end # function redfield

"Calculate the relative molar mass of POM in g/mol."
function rmm_pom(RC::Float64, RN::Float64, RP::Float64)
    rmm_CH2O = 30.031 # g/mol
    rmm_NH3 = 17.031 # g/mol
    rmm_H3PO4 = 97.994 # g/mol
    return RC*rmm_CH2O + RN*rmm_NH3 + RP*rmm_H3PO4
end # function rmm_pom

"Calculate the depth-dependent porosity."
function phi(phi0::Float64, phiInf::Float64, beta::Float64,
        depths::Array{Float64})
    return @. (phi0 - phiInf)*exp(-beta*depths) + phiInf
end # function phi

"Calculate the depth-dependent porosity - solid volume fraction."
phiS(phi::Array{Float64}) = 1.0 .- phi

"Calculate the tortuousity squared following Boudreau (1996, GCA)."
tort2(phi::Array{Float64}) = @. 1.0 - 2.0log(phi)

"Calculate the 1st derivative of the depth-dependent porosity w.r.t. depth."
function delta_phi(phi0::Float64, phiInf::Float64, beta::Float64,
        depths::Array{Float64})
    return @. -beta*(phi0 - phiInf)*exp(-beta*depths)
end # function delta_phi

"""Calculate the 1st derivative of the depth-dependent porosity - solid volume
fraction - w.r.t. depth.
"""
delta_phiS(delta_phi::Array{Float64}) = -delta_phi

"""Calculate the 1st derivative of the inverse of the tortuousity squared
w.r.t. depth.
"""
function delta_tort2i(delta_phi::Array{Float64}, phi::Array{Float64},
        tort2::Array{Float64})
    return @. 2.0delta_phi/(phi*tort2^2)
end # function delta_tort2i

"Calculate the diffusion-by-bioturbation constant coefficient."
D_bio_0(Fpoc::Float64) = @. 0.0232e-4*(1e2Fpoc)^0.85

"Calculate diffusion by bioturbation vs depth."
function D_bio(depths::Array{Float64}, D_bio_0::Float64, lambda_b::Float64,
        dO2_w::Float64)
    return @. D_bio_0*exp(-(depths/lambda_b)^2)*dO2_w/(dO2_w + 0.02)
end # function D_bio

"Calculate 1st derivative of diffusion by bioturbation w.r.t. depth."
function delta_D_bio(depths::Array{Float64}, D_bio::Array{Float64},
        lambda_b::Float64)
    return @. -2.0depths*D_bio/lambda_b^2
end # function delta_D_bio

"Calculate POC degradation parameter."
function krefractory(depths::Array{Float64,1}, D_bio_0::Float64)
    return @. 80.25D_bio_0*exp(-depths)
end # function krefractory

"Calculate bulk burial velocity at the sediment-water interface in m/a."
x0(Fpom::Float64, rho_pom::Float64, phiS_2::Float64) = Fpom/(rho_pom*phiS_2)

"Calculate bulk burial velocity at infinite depth in the sediment in m/a."
xinf(x0::Float64, phiS_2::Float64, phiS_e2::Float64) = x0*phiS_2/phiS_e2

"Calculate porewater burial velocity in m/a."
u(xinf::Float64, phi::Array{Float64}) = xinf*phi[end-1]./phi

"Calculate solid burial velocity in m/a."
w(xinf::Float64, phiS::Array{Float64}) = xinf*phiS[end-1]./phiS

"Calculate half of the cell Peclet number (Boudreau 1996 GCA, eq. 97)."
function Peh(w::Array{Float64}, z_res::Float64, D_bio::Array{Float64})
    return @. w*z_res/2.0D_bio
end # function Peh

"Calculate sigma (Boudreau 1996 GCA, eq. 96)."
sigma(Peh::Array{Float64}) = @. 1.0/tanh(Peh) - 1.0/(Peh)
function sigma(w::Array{Float64}, z_res::Float64, D_bio::Array{Float64})
    return sigma(Peh(w, z_res, D_bio))
end # function sigma

"Calculate T-dependent 'free solution' diffusion coeff. for O2 in m^2/a."
D_dO2(T::Float64) = 0.034862 + 0.001409T

"""Calculate T-dependent 'free solution' diffusion coeff. for tCO2 in m^2/a.
approximted to bicarbonate diffusion coefficient from Hulse et al (2018).
"""
D_dtCO2(T::Float64) = 0.015169 + 0.000793T

"Calculate alpha_0 parameter for irrigation (Archer et al. 2002)."
function alpha_0(Fpoc::Float64, dO2_w::Float64)
    return @. 11.0*(atan((1e2Fpoc*5.0 - 400.0)/400.0)/pi + 0.5) - 0.9 +
        20.0*(dO2_w/(dO2_w + 0.01))*exp(-dO2_w/0.01)*1e2Fpoc/(1e2Fpoc + 30.0)
end # function alpha_0

"Calculate alpha parameter for irrigation (Archer et al. 2002)."
function alpha(alpha_0::Float64, depths::Array{Float64}, lambda_i::Float64)
    return @. alpha_0*exp(-(depths/lambda_i)^2)
end # function alpha

"Calculate alpha parameter for irrigation (Archer et al. 2002)."
function alpha(Fpoc::Float64, dO2_w::Float64, depths::Array{Float64},
        lambda_i::Float64)
    return @. alpha_0(Fpoc, dO2_w)*exp(-(depths/lambda_i)^2)
end # function alpha

"Calculate APPW convenience term."
function APPW(w::Array{Float64}, delta_D_bio::Array{Float64},
        delta_phiS::Array{Float64}, D_bio::Array{Float64}, phiS::Array{Float64})
    return @. w - delta_D_bio - delta_phiS*D_bio/phiS
end # function APPW

"Calculate TR convenience term."
TR(z_res::Float64, tort2_2::Float64, dbl::Float64) = 2.0z_res*tort2_2/dbl

end # module Params
