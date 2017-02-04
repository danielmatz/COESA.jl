using COESA
using Base.Test

M0 = COESA.M0

# altitude conversions
@test isapprox(round(COESA.geopotential_altitude(0)), 0)
@test isapprox(round(COESA.geopotential_altitude(-5000)), -5004)
@test isapprox(round(COESA.geopotential_altitude(85_500)), 84_365)
@test isapprox(round(COESA.geopotential_altitude(1_000_000)), 864_071)

# out of bounds errors
@test_throws ErrorException atmosphere(prevfloat(-5_000.0))
@test_throws ErrorException atmosphere(nextfloat(1_000_000.0))

# b = 0, Z < 0
Z = -5_000

atmos = atmosphere(Z)
@test isapprox(signif(mean_molecular_weight(atmos), 5), 28.964)
@test isapprox(signif(temperature(atmos), 6), 320.676)
@test isapprox(signif(pressure(atmos), 5), 1.7776e5)
@test isapprox(signif(density(atmos), 5), 1.9311)
@test isapprox(signif(speed_of_sound(atmos), 5), 358.99)

# b = 0, Z = 0
Z = 0

atmos = atmosphere(Z)
@test isapprox(signif(mean_molecular_weight(atmos), 5), 28.964)
@test isapprox(signif(temperature(atmos), 6), 288.150)
@test isapprox(signif(pressure(atmos), 6), 101325)
@test isapprox(signif(density(atmos), 5), 1.2250)
@test isapprox(signif(speed_of_sound(atmos), 5), 340.29)

# b = 0, Z > 0
Z = 5_000

atmos = atmosphere(Z)
@test isapprox(signif(mean_molecular_weight(atmos), 5), 28.964)
@test isapprox(signif(temperature(atmos), 6), 255.676)
@test isapprox(signif(pressure(atmos), 5), 5.4048e4)
@test isapprox(signif(density(atmos), 5), 7.3643e-1)
@test isapprox(signif(speed_of_sound(atmos), 5), 320.55)

# b = 1
Z = 15_000

atmos = atmosphere(Z)
@test isapprox(signif(mean_molecular_weight(atmos), 5), 28.964)
@test isapprox(signif(temperature(atmos), 6), 216.650)
@test isapprox(signif(pressure(atmos), 4), signif(1.2111e4, 4))
@test isapprox(signif(density(atmos), 5), 1.9476e-1)
@test isapprox(signif(speed_of_sound(atmos), 5), 295.07)

# b = 2
Z = 25_000

atmos = atmosphere(Z)
@test isapprox(signif(mean_molecular_weight(atmos), 5), 28.964)
@test isapprox(signif(temperature(atmos), 6), 221.552)
@test isapprox(signif(pressure(atmos), 5), 2.5492e3)
@test isapprox(signif(density(atmos), 5), 4.0084e-2)
@test isapprox(signif(speed_of_sound(atmos), 5), 298.39)

# b = 3
Z = 40_000

atmos = atmosphere(Z)
@test isapprox(signif(mean_molecular_weight(atmos), 5), 28.964)
@test isapprox(signif(temperature(atmos), 6), 250.350)
@test isapprox(signif(pressure(atmos), 5), 2.8714e2)
@test isapprox(signif(density(atmos), 5), 3.9957e-3)
@test isapprox(signif(speed_of_sound(atmos), 5), 317.19)

# b = 4
Z = 50_000

atmos = atmosphere(Z)
@test isapprox(signif(mean_molecular_weight(atmos), 5), 28.964)
@test isapprox(signif(temperature(atmos), 6), 270.650)
@test isapprox(signif(pressure(atmos), 5), 7.9779e1)
@test isapprox(signif(density(atmos), 5), 1.0269e-3)
@test isapprox(signif(speed_of_sound(atmos), 5), 329.80)

# b = 5
Z = 60_000

atmos = atmosphere(Z)
@test isapprox(signif(mean_molecular_weight(atmos), 5), 28.964)
@test isapprox(signif(temperature(atmos), 6), 247.021)
@test isapprox(signif(pressure(atmos), 4), signif(2.1958e1, 4))
@test isapprox(signif(density(atmos), 5), 3.0968e-4)
@test isapprox(signif(speed_of_sound(atmos), 5), 315.07)

# b = 6
Z = 75_000

atmos = atmosphere(Z)
@test isapprox(signif(mean_molecular_weight(atmos), 5), 28.964)
@test isapprox(signif(temperature(atmos), 6), 208.399)
@test isapprox(signif(pressure(atmos), 5), 2.3881)
@test isapprox(signif(density(atmos), 5), 3.9921e-5)
@test isapprox(signif(speed_of_sound(atmos), 5), 289.40)

# b = 7
# The tables show values of M that have a discontinuity at Z = 86km, but our
# routine applies the suggested blending
Z = 85_000

atmos = atmosphere(Z)
@test isapprox(signif(mean_molecular_weight(atmos), 6), signif(M0 * 0.999694, 6))
@test isapprox(signif(temperature(atmos), 6), signif(188.893 / M0 * mean_molecular_weight(atmos), 6))
@test isapprox(signif(pressure(atmos), 5), 4.4568e-1)
@test isapprox(signif(density(atmos), 4), signif(8.2196e-6, 4))
@test isapprox(signif(speed_of_sound(atmos), 5), 275.52)

# b = 8, on a breakpoint
Z = 86_000

atmos = atmosphere(Z)
@test isapprox(signif(temperature(atmos), 5), 186.87)
@test isapprox(signif(pressure(atmos), 5), 3.7338e-1)
@test isapprox(signif(mean_molecular_weight(atmos), 4), 28.95)
@test isapprox(signif(density(atmos), 3), signif(6.958e-6, 3))
@test isapprox(signif(speed_of_sound(atmos), 5), 274.10)

# b = 8, in between breakpoints
Z = 86_500

atmos = atmosphere(Z)
@test isapprox(signif(temperature(atmos), 5), 186.87)
@test isapprox(signif(pressure(atmos), 5), 3.4163e-1)
@test isapprox(signif(mean_molecular_weight(atmos), 4), 28.95)
@test isapprox(signif(density(atmos), 4), 6.366e-6)
@test isapprox(signif(speed_of_sound(atmos), 5), 274.10)

# b = 8, in between breakpoints
Z = 100_000

atmos = atmosphere(Z)
@test isapprox(signif(temperature(atmos), 5), 195.08)
@test isapprox(signif(pressure(atmos), 2), signif(3.2011e-2, 2))
@test isapprox(signif(mean_molecular_weight(atmos), 4), 28.40)
@test isapprox(signif(density(atmos), 2), signif(5.604e-7, 2))
@test isapprox(signif(speed_of_sound(atmos), 5), 274.10)

# b = 9
Z = 115_000

atmos = atmosphere(Z)
@test isapprox(signif(temperature(atmos), 5), 300.00)
@test isapprox(signif(pressure(atmos), 5), 4.0096e-3)
@test isapprox(signif(mean_molecular_weight(atmos), 4), 26.68)
@test isapprox(signif(density(atmos), 4), 4.289e-8)
@test isapprox(signif(speed_of_sound(atmos), 5), 274.10)

# b = 10
Z = 200_000

atmos = atmosphere(Z)
@test isapprox(signif(temperature(atmos), 5), 854.56)
@test isapprox(signif(pressure(atmos), 5), 8.4736e-5)
@test isapprox(signif(mean_molecular_weight(atmos), 4), 21.30)
@test isapprox(signif(density(atmos), 3), signif(2.541e-10, 3))
@test isapprox(signif(speed_of_sound(atmos), 5), 274.10)

# b = 11
Z = 750_000

atmos = atmosphere(Z)
@test isapprox(signif(temperature(atmos), 5), 999.99)
@test isapprox(signif(pressure(atmos), 5), 2.2599e-8)
@test isapprox(signif(mean_molecular_weight(atmos), 3), 6.58)
@test isapprox(signif(density(atmos), 3), signif(1.788e-14, 3))
@test isapprox(signif(speed_of_sound(atmos), 5), 274.10)

# b = 11, in between breakpoints
Z = 985_000

atmos = atmosphere(Z)
@test isapprox(signif(temperature(atmos), 6), 1000.00)
@test isapprox(signif(pressure(atmos), 3), signif(7.9185e-9, 3))
@test isapprox(signif(mean_molecular_weight(atmos), 3), 3.99)
@test isapprox(signif(density(atmos), 3), signif(3.797e-15, 3))
@test isapprox(signif(speed_of_sound(atmos), 5), 274.10)

# b = 12
Z = 1_000_000

atmos = atmosphere(Z)
@test isapprox(signif(temperature(atmos), 6), 1000.00)
@test isapprox(signif(pressure(atmos), 5), 7.5138e-9)
@test isapprox(signif(mean_molecular_weight(atmos), 3), 3.94)
@test isapprox(signif(density(atmos), 4), 3.561e-15)
@test isapprox(signif(speed_of_sound(atmos), 5), 274.10)
