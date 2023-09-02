using COESA
using Test
using Aqua
using JET
using Unitful

quantity_type(T, units) = Quantity{T, dimension(units), typeof(units)}

M0 = COESA.M0

@testset "COESA" begin

    @testset "Altitude conversions" begin
        @test isapprox(round(COESA.geopotential_altitude(0)), 0)
        @test isapprox(round(COESA.geopotential_altitude(-5000)), -5004)
        @test isapprox(round(COESA.geopotential_altitude(85_500)), 84_365)
        @test isapprox(round(COESA.geopotential_altitude(1_000_000)), 864_071)
    end

    @testset "Altitude bounds" begin
        @test_throws ErrorException atmosphere(prevfloat(-5_000.0))
        @test_throws ErrorException atmosphere(nextfloat(1_000_000.0))
    end

    @testset "Z = -5,000" begin
        # b = 0, Z < 0
        atmos = atmosphere(-5_000)
        @test round(mean_molecular_weight(atmos), sigdigits = 5) == 28.964
        @test round(temperature(atmos), sigdigits = 6) == 320.676
        @test round(pressure(atmos), sigdigits = 5) == 1.7776e5
        @test round(density(atmos), sigdigits = 5) == 1.9311
        @test round(speed_of_sound(atmos), sigdigits = 5) == 358.99
        @test round(dynamic_viscosity(atmos), sigdigits = 5) == 1.9422e-5
    end

    @testset "Z = 0" begin
        # b = 0, Z = 0
        atmos = atmosphere(0)
        @test round(mean_molecular_weight(atmos), sigdigits = 5) == 28.964
        @test round(temperature(atmos), sigdigits = 6) == 288.150
        @test round(pressure(atmos), sigdigits = 6) == 101325
        @test round(density(atmos), sigdigits = 5) == 1.2250
        @test round(speed_of_sound(atmos), sigdigits = 5) == 340.29
        @test round(dynamic_viscosity(atmos), sigdigits = 5) == 1.7894e-5
    end

    @testset "Z = 5,000" begin
        # b = 0, Z > 0
        atmos = atmosphere(5_000)
        @test round(mean_molecular_weight(atmos), sigdigits = 5) == 28.964
        @test round(temperature(atmos), sigdigits = 6) == 255.676
        @test round(pressure(atmos), sigdigits = 5) == 5.4048e4
        @test round(density(atmos), sigdigits = 5) == 7.3643e-1
        @test round(speed_of_sound(atmos), sigdigits = 5) == 320.55
        @test round(dynamic_viscosity(atmos), sigdigits = 5) == 1.6282e-5
    end

    @testset "Z = 15,000" begin
        # b = 1
        atmos = atmosphere(15_000)
        @test round(mean_molecular_weight(atmos), sigdigits = 5) == 28.964
        @test round(temperature(atmos), sigdigits = 6) == 216.650
        @test round(pressure(atmos), sigdigits = 4) == round(1.2111e4, sigdigits = 4)
        @test round(density(atmos), sigdigits = 5) == 1.9476e-1
        @test round(speed_of_sound(atmos), sigdigits = 5) == 295.07
        @test round(dynamic_viscosity(atmos), sigdigits = 5) == 1.4216e-5
    end

    @testset "Z = 25,000" begin
        # b = 2
        atmos = atmosphere(25_000)
        @test round(mean_molecular_weight(atmos), sigdigits = 5) == 28.964
        @test round(temperature(atmos), sigdigits = 6) == 221.552
        @test round(pressure(atmos), sigdigits = 5) == 2.5492e3
        @test round(density(atmos), sigdigits = 5) == 4.0084e-2
        @test round(speed_of_sound(atmos), sigdigits = 5) == 298.39
        @test round(dynamic_viscosity(atmos), sigdigits = 5) == 1.4484e-5
    end

    @testset "Z = 40,000" begin
        # b = 3
        atmos = atmosphere(40_000)
        @test round(mean_molecular_weight(atmos), sigdigits = 5) == 28.964
        @test round(temperature(atmos), sigdigits = 6) == 250.350
        @test round(pressure(atmos), sigdigits = 5) == 2.8714e2
        @test round(density(atmos), sigdigits = 5) == 3.9957e-3
        @test round(speed_of_sound(atmos), sigdigits = 5) == 317.19
        @test round(dynamic_viscosity(atmos), sigdigits = 5) == 1.6009e-5
    end

    @testset "Z = 50,000" begin
        # b = 4
        atmos = atmosphere(50_000)
        @test round(mean_molecular_weight(atmos), sigdigits = 5) == 28.964
        @test round(temperature(atmos), sigdigits = 6) == 270.650
        @test round(pressure(atmos), sigdigits = 5) == 7.9779e1
        @test round(density(atmos), sigdigits = 5) == 1.0269e-3
        @test round(speed_of_sound(atmos), sigdigits = 5) == 329.80
        @test round(dynamic_viscosity(atmos), sigdigits = 5) == 1.7037e-5
    end

    @testset "Z = 60,000" begin
        # b = 5
        atmos = atmosphere(60_000)
        @test round(mean_molecular_weight(atmos), sigdigits = 5) == 28.964
        @test round(temperature(atmos), sigdigits = 6) == 247.021
        @test round(pressure(atmos), sigdigits = 4) == round(2.1958e1, sigdigits = 4)
        @test round(density(atmos), sigdigits = 5) == 3.0968e-4
        @test round(speed_of_sound(atmos), sigdigits = 5) == 315.07
        @test round(dynamic_viscosity(atmos), sigdigits = 5) == 1.5837e-5
    end

    @testset "Z = 75,000" begin
        # b = 6
        atmos = atmosphere(75_000)
        @test round(mean_molecular_weight(atmos), sigdigits = 5) == 28.964
        @test round(temperature(atmos), sigdigits = 6) == 208.399
        @test round(pressure(atmos), sigdigits = 5) == 2.3881
        @test round(density(atmos), sigdigits = 5) == 3.9921e-5
        @test round(speed_of_sound(atmos), sigdigits = 5) == 289.40
        @test round(dynamic_viscosity(atmos), sigdigits = 5) == 1.3759e-5
    end

    @testset "Z = 85,000" begin
        # b = 7
        # The tables show values of M that have a discontinuity at Z = 86km, but our
        # routine applies the suggested blending
        atmos = atmosphere(85_000)
        @test round(mean_molecular_weight(atmos), sigdigits = 6) == let
            M = M0
            Mscale = 0.999694
            round(M * Mscale, sigdigits = 6)
        end
        @test round(temperature(atmos), sigdigits = 6) == let
            T = 188.893
            Tscale = mean_molecular_weight(atmos) / M0
            round(T * Tscale, sigdigits = 6)
        end
        @test round(pressure(atmos), sigdigits = 5) == 4.4568e-1
        @test round(density(atmos), sigdigits = 4) == let
            ρ = 8.2196e-6
            M = M0
            Mscale = 0.999694
            T = 188.893
            Tscale = mean_molecular_weight(atmos) / M0
            round(ρ * Mscale / Tscale, sigdigits = 4)
        end
        @test round(speed_of_sound(atmos), sigdigits = 5) == 275.52
        @test round(dynamic_viscosity(atmos), sigdigits = 5) == let
            μ = 1.2647e-5
            T = 188.893
            S = 110.4
            Tscale = mean_molecular_weight(atmos) / M0
            round(μ * (T + S) * Tscale^(3/2) / (T * Tscale + S), sigdigits = 5)
        end
    end

    @testset "Z = 86,000" begin
        # b = 8, on a breakpoint
        atmos = atmosphere(86_000)
        @test round(temperature(atmos), sigdigits = 5) == 186.87
        @test round(pressure(atmos), sigdigits = 5) == 3.7338e-1
        @test round(mean_molecular_weight(atmos), sigdigits = 4) == 28.95
        @test round(density(atmos), sigdigits = 3) == round(6.958e-6, sigdigits = 3)
        @test round(speed_of_sound(atmos), sigdigits = 5) == 274.10
    end

    @testset "Z = 86,500" begin
        # b = 8, in between breakpoints
        atmos = atmosphere(86_500)
        @test round(temperature(atmos), sigdigits = 5) == 186.87
        @test round(pressure(atmos), sigdigits = 5) == 3.4163e-1
        @test round(mean_molecular_weight(atmos), sigdigits = 4) == 28.95
        @test round(density(atmos), sigdigits = 4) == 6.366e-6
        @test round(speed_of_sound(atmos), sigdigits = 5) == 274.10
        @test_throws ErrorException dynamic_viscosity(atmos)
    end

    @testset "Z = 100,000" begin
        # b = 8, in between breakpoints
        atmos = atmosphere(100_000)
        @test round(temperature(atmos), sigdigits = 5) == 195.08
        @test round(pressure(atmos), sigdigits = 2) == round(3.2011e-2, sigdigits = 2)
        @test round(mean_molecular_weight(atmos), sigdigits = 4) == 28.40
        @test round(density(atmos), sigdigits = 2) == round(5.604e-7, sigdigits = 2)
        @test round(speed_of_sound(atmos), sigdigits = 5) == 274.10
        @test_throws ErrorException dynamic_viscosity(atmos)
    end

    @testset "Z = 115,000" begin
        # b = 9
        atmos = atmosphere(115_000)
        @test round(temperature(atmos), sigdigits = 5) == 300.00
        @test round(pressure(atmos), sigdigits = 5) == 4.0096e-3
        @test round(mean_molecular_weight(atmos), sigdigits = 4) == 26.68
        @test round(density(atmos), sigdigits = 4) == 4.289e-8
        @test round(speed_of_sound(atmos), sigdigits = 5) == 274.10
        @test_throws ErrorException dynamic_viscosity(atmos)
    end

    @testset "Z = 200,000" begin
        # b = 10
        atmos = atmosphere(200_000)
        @test round(temperature(atmos), sigdigits = 5) == 854.56
        @test round(pressure(atmos), sigdigits = 5) == 8.4736e-5
        @test round(mean_molecular_weight(atmos), sigdigits = 4) == 21.30
        @test round(density(atmos), sigdigits = 3) == round(2.541e-10, sigdigits = 3)
        @test round(speed_of_sound(atmos), sigdigits = 5) == 274.10
        @test_throws ErrorException dynamic_viscosity(atmos)
    end

    @testset "Z = 750,000" begin
        # b = 11
        atmos = atmosphere(750_000)
        @test round(temperature(atmos), sigdigits = 5) == 999.99
        @test round(pressure(atmos), sigdigits = 5) == 2.2599e-8
        @test round(mean_molecular_weight(atmos), sigdigits = 3) == 6.58
        @test round(density(atmos), sigdigits = 3) == round(1.788e-14, sigdigits = 3)
        @test round(speed_of_sound(atmos), sigdigits = 5) == 274.10
        @test_throws ErrorException dynamic_viscosity(atmos)
    end

    @testset "Z = 985,000" begin
        # b = 11, in between breakpoints
        atmos = atmosphere(985_000)
        @test round(temperature(atmos), sigdigits = 6) == 1000.00
        @test round(pressure(atmos), sigdigits = 3) == round(7.9185e-9, sigdigits = 3)
        @test round(mean_molecular_weight(atmos), sigdigits = 3) == 3.99
        @test round(density(atmos), sigdigits = 3) == round(3.797e-15, sigdigits = 3)
        @test round(speed_of_sound(atmos), sigdigits = 5) == 274.10
        @test_throws ErrorException dynamic_viscosity(atmos)
    end

    @testset "Z = 1,000,000" begin
        # b = 12
        atmos = atmosphere(1_000_000)
        @test round(temperature(atmos), sigdigits = 6) == 1000.00
        @test round(pressure(atmos), sigdigits = 5) == 7.5138e-9
        @test round(mean_molecular_weight(atmos), sigdigits = 3) == 3.94
        @test round(density(atmos), sigdigits = 4) == 3.561e-15
        @test round(speed_of_sound(atmos), sigdigits = 5) == 274.10
        @test_throws ErrorException dynamic_viscosity(atmos)
    end

    @testset "Unitful Extension" begin
        atmos = atmosphere(0u"m")
        @test altitude(atmos) == 0u"m"
        @test round(quantity_type(Float64, u"kg/kmol"), mean_molecular_weight(atmos), sigdigits = 5) == 28.964u"kg/kmol"
        @test round(quantity_type(Float64, u"K"), temperature(atmos), sigdigits = 6) == 288.150u"K"
        @test round(quantity_type(Float64, u"Pa"), pressure(atmos), sigdigits = 6) == 101325u"Pa"
        @test round(quantity_type(Float64, u"kg/m^3"), density(atmos), sigdigits = 5) == 1.2250u"kg/m^3"
        @test round(quantity_type(Float64, u"m/s"), speed_of_sound(atmos), sigdigits = 5) == 340.29u"m/s"
        @test round(quantity_type(Float64, u"N*s/m^2"), dynamic_viscosity(atmos), sigdigits = 5) == 1.7894e-5u"N*s/m^2"
    end

    @testset "Code Quality" begin
        @testset "Aqua" begin
            Aqua.test_all(COESA, project_extras = false)
        end
        @testset "JET" begin
            test_package("COESA")
        end
    end

end
