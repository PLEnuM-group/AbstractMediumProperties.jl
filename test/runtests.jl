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
    @test scattering_function(wavelength, medium) isa Number
end

struct MockMediumProperties <: MediumProperties end

AbstractMediumProperties.pressure(medium::MockMediumProperties) = 1.
AbstractMediumProperties.temperature(medium::MockMediumProperties) = 2.
AbstractMediumProperties.material_density(medium::MockMediumProperties) = 3.
AbstractMediumProperties.mean_scattering_angle(medium::MockMediumProperties) = 4.
AbstractMediumProperties.scattering_function(medium::MockMediumProperties) = π
AbstractMediumProperties.scattering_length(wavelength, medium::MockMediumProperties) = 5.
AbstractMediumProperties.absorption_length(wavelength, medium::MockMediumProperties) = 6.
AbstractMediumProperties.group_refractive_index(wavelength, medium::MockMediumProperties) = 1.2
AbstractMediumProperties.phase_refractive_index(wavelength, medium::MockMediumProperties) = 1.1
AbstractMediumProperties.dispersion(wavelength, medium::MockMediumProperties) = 0.001


@testset "AbstractMediumProperties.jl" begin
    medium = MockMediumProperties()
    test_interface(medium)

    wavelength = 400.
    @test group_velocity(400., medium) ≈ 0.45423099697
    @test cherenkov_angle(400., medium) ≈ 0.429699666
end
