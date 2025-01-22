using AbstractMediumProperties
using Test

function test_interface(medium::MediumProperties)
    wavelength = 400.

    @test pressure(medium) isa Number
    @test temperature(medium) isa Number
    @test material_density(medium) isa Number
    @test mean_scattering_angle(medium) isa Number
    @test scattering_length(wavelength, medium) isa Number
    @test absorption_length(wavelength, medium) isa Number
    @test group_refractive_index(wavelength, medium) isa Number
    @test phase_refractive_index(wavelength, medium) isa Number
    @test dispersion(wavelength, medium) isa Number
    @test cherenkov_angle(wavelength, medium) isa Number
    @test group_velocity(wavelength, medium) isa Number
    @test scattering_function(medium) isa Number
    @test radiation_length(medium) isa Number
end

# Define a concrete MediumProperties subtype for testing
struct MockMediumProperties <: MediumProperties
    salinity::Float64
    temperature::Float64
    pressure::Float64
end

AbstractMediumProperties.pressure(medium::MockMediumProperties) = 1.
AbstractMediumProperties.temperature(medium::MockMediumProperties) = 2.
AbstractMediumProperties.material_density(medium::MockMediumProperties) = 3.
AbstractMediumProperties.mean_scattering_angle(medium::MockMediumProperties) = 4.
AbstractMediumProperties.scattering_function(medium::MockMediumProperties) = π
AbstractMediumProperties.scattering_length(wavelength, medium::MockMediumProperties) = 5.
AbstractMediumProperties.absorption_length(wavelength, medium::MockMediumProperties) = 6.
AbstractMediumProperties.group_refractive_index(wavelength, medium::MockMediumProperties) = 1.2
AbstractMediumProperties.radiation_length(medium::MockMediumProperties) = 30

function AbstractMediumProperties.phase_refractive_index(wavelength, medium::MockMediumProperties)
    return refractive_index_fry(wavelength; salinity=medium.salinity, temperature=medium.temperature, pressure=medium.pressure)
end

function AbstractMediumProperties.dispersion(wavelength, medium::MockMediumProperties)
    return dispersion_fry(wavelength; salinity=medium.salinity, temperature=medium.temperature, pressure=medium.pressure)
end


@testset "AbstractMediumProperties.jl" begin
    medium = MockMediumProperties(35.0, 15.0, 1.)
    test_interface(medium)
end

# Test phase_refractive_index function
@testset "Phase Refractive Index" begin
    medium = MockMediumProperties(35.0, 15.0, 1.0)
    wavelength = 500.0 # nm
    n_phase = phase_refractive_index(wavelength, medium)
    @test isapprox(n_phase, 1.34, atol=0.01)
end

# Test dispersion function
@testset "Dispersion" begin
    medium = MockMediumProperties(35.0, 15.0, 1.0)
    wavelength = 500.0 # nm
    disp = dispersion(wavelength, medium)
    @test isapprox(disp, -0.0001, atol=0.0001)
end

# Test group_refractive_index function
@testset "Group Refractive Index" begin
    medium = MockMediumProperties(35.0, 15.0, 1.0)
    wavelength = 500.0 # nm
    n_group = group_refractive_index(wavelength, medium)
    @test isapprox(n_group, 1.2, atol=0.01)
end

# Test group_velocity function
@testset "Group Velocity" begin
    medium = MockMediumProperties(35.0, 15.0, 1.0)
    wavelength = 500.0 # nm
    v_group = group_velocity(wavelength, medium)
    @test isapprox(v_group, 0.22, atol=0.01)
end

# Test cherenkov_angle function
@testset "Cherenkov Angle" begin
    medium = MockMediumProperties(35.0, 15.0, 1.0)
    wavelength = 500.0 # nm
    angle = cherenkov_angle(wavelength, medium)
    @test isapprox(angle, 0.73, atol=0.01)
end
