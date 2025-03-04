export DIPPR105
export refractive_index_fry
export dispersion_fry
export sca_len_part_conc


# Scattering functions
export hg_scattering_func
export sl_scattering_func
export generic_scattering_function
export mixed_hg_sl_scattering_func
export mixed_hg_sl_scattering_func_ppc
export mixed_hg_generic_scattering_func

using Polynomials


"""
DIPPR105Params

Parameters for the DIPPR105 formula
"""
@kwdef struct DIPPR105Params
A::Float64
B::Float64
C::Float64
D::Float64
end

# DIPR105 Parameters from DDB
const DDBDIPR105Params = DIPPR105Params(A=0.14395, B=0.0112, C=649.727, D=0.05107)

"""
DIPPR105(temperature::Real, params::DIPPR105Params=DDBDIPR105Params)

Use DPPIR105 formula to calculate water density as function of temperature.
temperature in K.

Reference: http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe?component=Water

Returns density in kg/m^3
"""
function DIPPR105(temperature::Real, params::DIPPR105Params=DDBDIPR105Params)
    return params.A / (params.B^(1 + (1 - temperature / params.C)^params.D))
end



