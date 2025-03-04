using CherenkovMediumBase
using Test
using Random

Random.seed!(1234)


function test_interface(medium::MediumProperties)
    wavelength = 400.

    @test pressure(medium) isa Number
    @test temperature(medium) isa Number
    @test material_density(medium) isa Number
    @test mean_scattering_angle(medium) isa Number
    @test scattering_length(medium, wavelength) isa Number
    @test absorption_length(medium, wavelength) isa Number
    @test group_refractive_index(medium, wavelength) isa Number
    @test phase_refractive_index(medium, wavelength) isa Number
    @test dispersion(medium, wavelength) isa Number
    @test cherenkov_angle(medium, wavelength) isa Number
    @test group_velocity(medium, wavelength) isa Number
    @test sample_scattering_function(medium) isa Number

end

# Define a dispersion model



# Define a concrete MediumProperties subtype for testing
struct MockMediumProperties <: MediumProperties
    salinity::Float64
    temperature::Float64
    pressure::Float64
    dispersion_model::QuanFryDispersion
    scattering_model::KopelevichScatteringModel
    MockMediumProperties(salinity, temperature, pressure) = new(salinity, temperature, pressure, QuanFryDispersion(salinity, temperature, pressure), KopelevichScatteringModel(HenyeyGreenStein(0.8), 0.1, 0.2))
end


CherenkovMediumBase.pressure(medium::MockMediumProperties) = medium.pressure
CherenkovMediumBase.temperature(medium::MockMediumProperties) = medium.temperature
CherenkovMediumBase.material_density(medium::MockMediumProperties) = 3.
CherenkovMediumBase.mean_scattering_angle(medium::MockMediumProperties) = 4.
CherenkovMediumBase.absorption_length(medium::MockMediumProperties, wavelength) = 6.
CherenkovMediumBase.sample_scattering_function(medium::MockMediumProperties) =  rand(medium.scattering_model.scattering_function)


function CherenkovMediumBase.phase_refractive_index(medium::MockMediumProperties, wavelength)
    return phase_refractive_index(medium.dispersion_model, wavelength)
end

function CherenkovMediumBase.dispersion(medium::MockMediumProperties, wavelength)
    return dispersion(medium.dispersion_model, wavelength)
end

function CherenkovMediumBase.scattering_length(medium::MockMediumProperties, wavelength)
    return scattering_length(medium.scattering_model, wavelength)
end


@testset "CherenkovMediumBase.jl" begin
    medium = MockMediumProperties(35.0, 15.0, 1.)
    test_interface(medium)
end

@testset "Dispersion Models" begin
    @testset "QuanFry Dispersion" begin
        disp_model = QuanFryDispersion(35.0, 15.0, 1.0)

        @testset "Phase Refractive Index" begin
            wavelength = 500.0 # nm
            n_phase = phase_refractive_index(disp_model, wavelength)
            @test isapprox(n_phase, 1.34, atol=0.01)
        end

        @testset "Dispersion" begin
            wavelength = 500.0 # nm
            disp = dispersion(disp_model, wavelength)
            @test isapprox(disp, -0.0001, atol=0.0001)
        end
        
        # Test group_refractive_index function
        @testset "Group Refractive Index" begin
            medium = MockMediumProperties(35.0, 15.0, 1.0)
            wavelength = 500.0 # nm
            n_group = group_refractive_index(medium, wavelength)
            @test isapprox(n_group, 1.368, atol=0.01)
        end

        # Test group_velocity function
        @testset "Group Velocity" begin
            medium = MockMediumProperties(35.0, 15.0, 1.0)
            wavelength = 500.0 # nm
            v_group = group_velocity(medium, wavelength)
            @test isapprox(v_group, 0.22, atol=0.01)
        end

        # Test cherenkov_angle function
        @testset "Cherenkov Angle" begin
            medium = MockMediumProperties(35.0, 15.0, 1.0)
            wavelength = 500.0 # nm
            angle = cherenkov_angle(medium, wavelength)
            @test isapprox(angle, 0.73, atol=0.01)
        end
    end
end

@testset "Scattering" begin

    @testset "Scattering Functions" begin
        @testset "HenyeyGreenStein" begin
            scattering_function = HenyeyGreenStein(0.95)
            @test rand(scattering_function) ≈ 0.99014 atol=1e-5
        end

        @testset "EinsteinSmoluchowsky" begin
            scattering_function = EinsteinSmoluchowsky(0.835)
            @test rand(scattering_function) ≈ -0.41788 atol=1e-5
        end    
    
    end

    @testset "KopelevichScatteringModel" begin
        scattering_model = KopelevichScatteringModel(HenyeyGreenStein(0.95), 7.5E-3, 7.5E-3)

        # Test scattering_length function
        @testset "scattering_length" begin
            @test scattering_length(scattering_model, 400.0) ≈ 37.69 atol=1e-2
        end
    end
end