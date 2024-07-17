module AbstractMediumProperties

export MediumProperties
export pressure, temperature, material_density, mean_scattering_angle
export scattering_length, absorption_length, group_refractive_index, phase_refractive_index
export dispersion, cherenkov_angle, group_velocity, radiation_length
export scattering_function



const c_vac_m_ns = 0.299792458

abstract type MediumProperties end

_not_implemented(medium) = error("Not implemented for type $(typeof(medium))")

"""
    pressure(::MediumProperties) 

This function returns the pressure for a given medium.
"""
pressure(medium::MediumProperties) = _not_implemented(medium)
temperature(medium::MediumProperties) = _not_implemented(medium) 
material_density(medium::MediumProperties) = _not_implemented(medium) 
mean_scattering_angle(medium::MediumProperties) = _not_implemented(medium)
radiation_length(medium::MediumProperties) = _not_implemented(medium)


"""
    scattering_function(medium::MediumProperties)

Return a scattering angle sampled from the scattering function of the medium.

"""
scattering_function(medium::MediumProperties) = _not_implemented(medium)

"""
    scattering_length(wavelength, medium::MediumProperties) 
Return scattering length at `wavelength` in units m.

`wavelength` is expected to be in units nm. Returned length is in units m.
"""
scattering_length(wavelength, medium::MediumProperties) = _not_implemented(medium)


"""
    absorption_length(wavelength, medium::MediumProperties)
Return absorption length at `wavelength` in units m.

`wavelength` is expected to be in units nm.
"""
absorption_length(wavelength, medium::MediumProperties) = _not_implemented(medium)





"""
    phase_refractive_index(wavelength, medium)
Return the phase refractive index at `wavelength`.

`wavelength` is expected to be in units nm.
"""
phase_refractive_index(wavelength, medium::MediumProperties) = _not_implemented(medium)


"""
    dispersion(wavelength, medium)
Return the dispersion dn/dλ at `wavelength` in units 1/nm.

`wavelength` is expected to be in units nm.
"""
dispersion(wavelength, medium::MediumProperties) = _not_implemented(medium)


"""
    cherenkov_angle(wavelength, medium::MediumProperties)
Calculate the cherenkov angle (in rad) for `wavelength`.

`wavelength` is expected to be in units nm.
"""
function cherenkov_angle(wavelength, medium::MediumProperties)
    return acos(one(typeof(wavelength)) / phase_refractive_index(wavelength, medium))
end

"""
    group_velocity(wavelength, medium)
Return the group_velocity in m/ns at `wavelength`.

`wavelength` is expected to be in units nm.
"""
function group_velocity(wavelength, medium::MediumProperties)
    global c_vac_m_ns
    T = typeof(wavelength)

    # Explicitely convert everything to the type of wavelength
    # This ise useful to avoud double precision if this function is called in a CUDA kernel

    ref_ix::T = phase_refractive_index(wavelength, medium)
    λ_0::T = ref_ix * wavelength
    T(c_vac_m_ns) / (ref_ix - λ_0 * dispersion(wavelength, medium))
end

"""
    group_refractive_index(wavelength, medium)
Return the group refractive index at `wavelength`.

`wavelength` is expected to be in units nm.
"""
function group_refractive_index(wavelength, medium)
    global c_vac_m_ns

    T = typeof(wavelength)

    return T(c_vac_m_ns)/group_velocity(wavelength, medium)
end

include("utils.jl")


end # module
