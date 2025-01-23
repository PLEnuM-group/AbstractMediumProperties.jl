import Distributions
using Distributions: Sampleable, Univariate, Continuous
using Random
using Polynomials: fit, Polynomial, ImmutablePolynomial

export AbstractScatteringFunction
export AbstractScatteringModel
export HenyeyGreenStein
export PolynomialScatteringFunction
export EinsteinSmoluchowsky
export SimplifiedLiu
export TwoComponentScatteringModel
export MixedHGES
export MixedHGSL

export KopelevichScatteringModel

abstract type AbstractScatteringFunction <: Sampleable{Univariate, Continuous} end

rand(rng::AbstractRNG, s::AbstractScatteringFunction) = _not_implemented(s)


struct HenyeyGreenStein{T} <: AbstractScatteringFunction
    g::T # Mean scattering angle
end

"""
    _hg_scattering_func(g::Real)

Henyey-Greenstein scattering in one plane.

# Arguments
- `g::Real`: mean scattering angle

# Returns
- `typeof(g)` cosine of a scattering angle sampled from the distribution

"""
function _hg_scattering_func(rng::AbstractRNG, g::T) where {T <: Real}
    eta = Base.rand(rng, T)
    costheta::T = (1 / (2 * g) * (1 + g^2 - ((1 - g^2) / (1 + g * (2 * eta - 1)))^2))
    return clamp(costheta, T(-1), T(1))
end

Base.rand(rng::AbstractRNG, s::HenyeyGreenStein) = _hg_scattering_func(rng, s.g)

struct SimplifiedLiu{T} <: AbstractScatteringFunction
    g::T # Mean scattering angle
end

"""
    sl_scattering_func(g::Real)
Simplified-Liu scattering angle function.
Implementation from: https://user-web.icecube.wisc.edu/~dima/work/WISC/ppc/spice/new/paper/a.pdf

# Arguments
- `g::Real`: mean scattering angle
"""
function sl_scattering_func(rng::AbstractRNG, g::T) where {T <: Real}
    eta = Base.rand(rng, T)
    beta = (1-g) / (1+g)
    costheta::T = 2 * eta^beta - 1
    return clamp(costheta, T(-1), T(1))
end

Base.rand(rng::AbstractRNG, s::SimplifiedLiu) = sl_scattering_func(rng, s.g)




struct PolynomialScatteringFunction{T, P <: ImmutablePolynomial{T}} <: AbstractScatteringFunction
    poly::P
end

function Base.rand(rng::AbstractRNG, s::PolynomialScatteringFunction{T}) where {T}
    eta = Base.rand(rng, T)
    return clamp(s.poly(eta), T(-1), T(1))
end

"""
    es_scattering(cos_theta::T, b::T) where {T<:Real}

Einstein-Smoluchowsky Scattering PDF.
"""
function es_scattering(cos_theta::T, b::T) where {T<:Real}
    a = 1/(4*pi) * 1 / (1+b/3)
    return a*(1+b*cos_theta^2)
end

"""
    es_scattering_integral(cos_theta::T, b::T) where {T<:Real}

Anti-derivative of the ES scattering function.
"""
function es_scattering_integral(cos_theta::T, b::T) where {T<:Real}
    a = 1/(4*pi) * 1 / (1+b/3)
    return a*cos_theta*(1 + (b*cos_theta^2)/3) * 2*pi
end

"""
    es_scattering_integral(cos_theta::T, b::T) where {T<:Real}

Integral of ES scattering function from -1 to cos_theta
"""
function es_scattering_cumulative(cos_theta::T, b::T) where {T<:Real}
    return es_scattering_integral(cos_theta, b) - es_scattering_integral(-1., b)
end

"""
    make_inverse_es_polynomial(b::T) where {T<:Real}

Make a 3rd order polynomial that approximates the inverse of the ES scattering function.
"""
function make_inverse_es_polynomial(b::T) where {T<:Real}
    cos_theta = -1:0.01:1
    es_poly = fit(Polynomial, es_scattering_cumulative.(cos_theta, b), cos_theta, 3)
    return ImmutablePolynomial(Tuple(T.(collect(es_poly))))
end

function EinsteinSmoluchowsky(b::T) where {T}
    poly = make_inverse_es_polynomial(b)
    return PolynomialScatteringFunction(poly)
end

struct TwoComponentScatteringModel{F1<:AbstractScatteringFunction, F2<:AbstractScatteringFunction} <: AbstractScatteringFunction
    f1::F1
    f2::F2
    fraction::Real
end


MixedHGES(g, b, fraction) = TwoComponentScatteringModel(HenyeyGreenStein(g), EinsteinSmoluchowsky(b), fraction)
MixedHGSL(g, fraction) = TwoComponentScatteringModel(HenyeyGreenStein(g), SimplifiedLiu(g), fraction)


function Base.rand(rng::AbstractRNG, s::TwoComponentScatteringModel)
    choice = Base.rand(rng, Float64)
    if choice < s.fraction
        return Base.rand(rng, s.f1)
    end
    return Base.rand(rng, s.f2)
end


"""
    sca_len_part_conc(wavelength; vol_conc_small_part, vol_conc_large_part)

Calculates the scattering length (in m) for a given wavelength based on concentrations of
small (`vol_conc_small_part`) and large (`vol_conc_large_part`) particles.
wavelength is given in nm, vol_conc_small_part and vol_conc_large_part in ppm

C.D. Mobley "Light and Water", ISBN 0-12-502750-8, pag. 119. 
"""
@inline function _sca_len_part_conc(
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


abstract type AbstractScatteringModel end
scattering_length(model::AbstractScatteringModel, wavelength::Real) = _not_implemented(model)

struct KopelevichScatteringModel{F<:AbstractScatteringFunction, T} <: AbstractScatteringModel
    scattering_function::F
    vol_conc_large_part::T
    vol_conc_small_part::T
end

function scattering_length(model::KopelevichScatteringModel, wavelength::Real)
    return _sca_len_part_conc(
        wavelength,
        vol_conc_small_part=model.vol_conc_small_part,
        vol_conc_large_part=model.vol_conc_large_part
    )
end



"""
    mixed_hg_es_scattering_func(g::Real, hg_fraction::Real, poly_coeffs::NTuple{4, Real})

Mixture model of HG and generic scattering function.
"""
function mixed_hg_generic_scattering_func(g::T, hg_fraction::T, poly_coeffs::NTuple{4, T}) where {T<:Real}
    choice = rand(T)
    if choice < hg_fraction
        return hg_scattering_func(g)
    end
    return generic_scattering_function(poly_coeffs...)
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

"""
    mixed_hg_sl_scattering_func_ppc(g, hg_fraction)
Mixture model of HG and SL as implemented in PPC
"""
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


