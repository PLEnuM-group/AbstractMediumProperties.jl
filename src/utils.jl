
export DIPPR105
export refractive_index_fry
export dispersion_fry
export sca_len_part_conc
export calc_quan_fry_params

export hg_scattering_func, sl_scattering_func, mixed_hg_sl_scattering_func, mixed_hg_sl_scattering_func_ppc

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

"""
    calc_quan_fry_params(salinity::Real, temperature::Real, pressure::Real)

Helper function to get the parameters for the Quan & Fry formula as function of
salinity, temperature and pressure.
"""
function calc_quan_fry_params(
    salinity::Real,
    temperature::Real,
    pressure::Real)

    n0 = 1.31405
    n1 = 1.45e-5
    n2 = 1.779e-4
    n3 = 1.05e-6
    n4 = 1.6e-8
    n5 = 2.02e-6
    n6 = 15.868
    n7 = 0.01155
    n8 = 0.00423
    n9 = 4382
    n10 = 1.1455e6

    a01 = (
        n0
        +
        (n2 - n3 * temperature + n4 * temperature^2) * salinity
        -
        n5 * temperature^2
        +
        n1 * pressure
    )
    a2 = n6 + n7 * salinity - n8 * temperature
    a3 = -n9
    a4 = n10

    return a01, a2, a3, a4
end

"""
    refractive_index_fry(wavelength, salinity, temperature, pressure)

The phase refractive index of sea water according to a model
from Quan & Fry.

wavelength is given in nm, salinity in permille, temperature in °C and pressure in atm

The original model is taken from:
X. Quan, E.S. Fry, Appl. Opt., 34, 18 (1995) 3477-3480.

An additional term describing pressure dependence was included according to:
Wolfgang H.W.A. Schuster, "Measurement of the Optical Properties of the Deep
Mediterranean - the ANTARES Detector Medium.",
PhD thesis (2002), St. Catherine's College, Oxford
downloaded Jan 2011 from: http://www.physics.ox.ac.uk/Users/schuster/thesis0098mmjhuyynh/thesis.ps

Adapted from clsim (©Claudio Kopper)
"""
function refractive_index_fry(
    wavelength::T;
    salinity::Real,
    temperature::Real,
    pressure::Real) where {T<:Real}
    refractive_index_fry(wavelength, T.(calc_quan_fry_params(salinity, temperature, pressure)))
end

function refractive_index_fry(
    wavelength::Real,
    quan_fry_params::Tuple{U,U,U,U}
) where {U<:Real}

    a01, a2, a3, a4 = quan_fry_params
    x = one(wavelength) / wavelength
    x2 = x * x
    # a01 + x*a2 + x^2 * a3 + x^3 * a4
    return oftype(wavelength, a01 + x*a2 + x2*a3 + x2 * x * a4)
end


"""
dispersion_fry(
    wavelength::T;
    salinity::Real,
    temperature::Real,
    pressure::Real) where {T <: Real}

    Calculate the dispersion (dn/dλ) for the Quan & Fry model.
    wavelength is given in nm, salinity in permille, temperature in °C and pressure in atm
"""
function dispersion_fry(
    wavelength::T;
    salinity::Real,
    temperature::Real,
    pressure::Real) where {T<:Real}
    dispersion_fry(wavelength, T.(calc_quan_fry_params(salinity, temperature, pressure)))
end

function dispersion_fry(wavelength::T, quan_fry_params::NTuple{4,<:Number}) where {T<:Number}
    a2, a3, a4 = quan_fry_params
    x = one(T) / wavelength

    return T(a2 + T(2) * x * a3 + T(3) * x^2 * a4) * T(-1) / wavelength^2
end




"""
    sca_len_part_conc(wavelength; vol_conc_small_part, vol_conc_large_part)

Calculates the scattering length (in m) for a given wavelength based on concentrations of
small (`vol_conc_small_part`) and large (`vol_conc_large_part`) particles.
wavelength is given in nm, vol_conc_small_part and vol_conc_large_part in ppm


Adapted from clsim ©Claudio Kopper
"""
@inline function sca_len_part_conc(
    wavelength::T;
    vol_conc_small_part::Real,
    vol_conc_large_part::Real) where {T<:Real}

    ref_wlen::T = 550  # nm
    x::T = ref_wlen / wavelength

    sca_coeff = (
        T(0.0017) * x^T(4.3)
        + T(1.34) * vol_conc_small_part * x^T(1.7)
        + T(0.312) * vol_conc_large_part * x^T(0.3)
    )

    return T(1 / sca_coeff)

end


"""
    hg_scattering_func(g::Real)

CUDA-optimized version of Henyey-Greenstein scattering in one plane.

# Arguments
- `g::Real`: mean scattering angle

# Returns
- `typeof(g)` cosine of a scattering angle sampled from the distribution

"""
@inline function hg_scattering_func(g::T) where {T <: Real}
    """Henyey-Greenstein scattering in one plane."""
    eta = rand(T)
    costheta::T = (1 / (2 * g) * (1 + g^2 - ((1 - g^2) / (1 + g * (2 * eta - 1)))^2))
    #costheta::T = (1 / (2 * g) * (fma(g, g, 1) - (fma(-g, g, 1) / (fma(g, (fma(2, eta, -1)), 1)))^2))
    return clamp(costheta, T(-1), T(1))



end

"""
    sl_scattering_func(g::Real)
Simplified-Liu scattering angle function.
Implementation from: https://user-web.icecube.wisc.edu/~dima/work/WISC/ppc/spice/new/paper/a.pdf

# Arguments
- `g::Real`: mean scattering angle
"""
function sl_scattering_func(g::T) where {T <: Real}
    eta = rand(T)
    beta = (1-g) / (1+g)
    costheta::T = 2 * eta^beta - 1
    return clamp(costheta, T(-1), T(1))
end



"""
    mixed_hg_sl_scattering_func(g, hg_fraction)
Mixture model of HG and SL.

# Arguments
- `g::Real`: mean scattering angle
- `hg_fraction::Real`: mixture weight of the HG component
"""
function mixed_hg_sl_scattering_func(g::Real, hg_fraction::Real)
    choice = rand()
    if choice < hg_fraction
        return hg_scattering_func(g)
    end
    return sl_scattering_func(g)
end

function mixed_hg_sl_scattering_func_ppc(g::T, hg_fraction::T) where {T <:Real}
    xi = rand()
    sf = hg_fraction
    gr::T = (1-g)/(1+g)
	if(xi>sf)
	  xi=(1-xi)/(1-sf)
	  xi=2*xi-1
	  if(g!=0)
	    ga::T=(1-g*g)/(1+g*xi)
	    xi=(1+g*g-ga*ga)/(2*g)
      end
	else
	  xi/=sf
	  xi=2*xi^gr-1
    end
    return clamp(xi, T(-1), T(1))
end

