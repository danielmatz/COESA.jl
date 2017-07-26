module COESA

export atmosphere, density, temperature, pressure, speed_of_sound, mean_molecular_weight

const r0 = 6356766.0 # (m), effective Earth radius at 45 deg N latitude
const g0 = 9.80665 # (m / s²) or (m² / s² m')
const M0 = 28.9644 # (kg / kmol)
const Rstar = 8.31432e3 # (N m / kmol K)
const γ = 1.4

geopotential_altitude(Z) = r0 * Z / (r0 + Z)

const Hb = [0, 11, 20, 32, 47, 51, 71] * 1000 # (m')
const Lmb = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0] / 1000 # (K / m')

function findb(H)
    i = 1
    while i < length(Hb) && H > Hb[i + 1]
        i += 1
    end
    b = i - 1
    b
end

const Tmb = let
    Tmb = [288.15] # (K)
    for i in 1:(length(Hb) - 1)
        push!(Tmb, Tmb[i] + Lmb[i] * (Hb[i + 1] - Hb[i]))
    end
    Tmb
end
function Tm(H)
    b = findb(H)
    i = b + 1
    Tmb[i] + Lmb[i] * (H - Hb[i])
end
temperature_lower(H, M) = Tm(H) / M0 * M

const Pb = let
    Pb = [101325.0]
    for i in 1:(length(Hb) - 1)
        if Lmb[i] == 0
            push!(Pb, Pb[i] * exp(-g0 * M0 * (Hb[i + 1] - Hb[i]) / (Rstar * Tmb[i])))
        else
            push!(Pb, Pb[i] * (Tmb[i] / (Tmb[i] + Lmb[i] * (Hb[i + 1] - Hb[i]))) ^ (g0 * M0 / (Rstar * Lmb[i])))
        end
    end
    Pb
end
function pressure_lower(H)
    b = findb(H)
    i = b + 1
    if Lmb[i] == 0
        return Pb[i] * exp(-g0 * M0 * (H - Hb[i]) / (Rstar * Tmb[i]))
    else
        return Pb[i] * (Tmb[i] / (Tmb[i] + Lmb[i] * (H - Hb[i]))) ^ (g0 * M0 / (Rstar * Lmb[i]))
    end
end

speed_of_sound_lower(T, M) = sqrt(γ * Rstar * T / M)

const Ztable = [80.0, 80.5, 81.0, 81.5, 82.0, 82.5, 83.0, 83.5, 84.0, 84.5, 85.0, 85.5, 86.0] * 1000 # (m)
const Mratiotable = [1.0, 0.999996, 0.999989, 0.999971, 0.999941, 0.999909, 0.999870, 0.999829, 0.999786, 0.999741, 0.999694, 0.999641, 0.999579]
function mean_molecular_weight_ratio_lower(Z)
    Z < Ztable[1] && return 1.0
    Z > Ztable[end] && error("altitude above maximum value in table")
    i = 1
    while i < length(Ztable) && Z > Ztable[i + 1]
        i += 1
    end
    Mratiotable[i] + (Mratiotable[i + 1] - Mratiotable[i]) / (Ztable[i + 1] - Ztable[i]) * (Z - Ztable[i])
end
mean_molecular_weight_lower(Z) = M0 * mean_molecular_weight_ratio_lower(Z)

function temperature_upper(Z)
    if Z <= 91_000
        return 186.8673 # (K)
    elseif Z <= 110_000
        Tc = 263.1905 # (K)
        A = -76.3232 # (K)
        a = -19.9429 * 1000 # (m)
        return Tc + A * sqrt(1 - ((Z - 91_000) / a) ^ 2)
    elseif Z <= 120_000
        T9 = 240 # (K)
        LK9 = 12 / 1000 # (K / m)
        Z9 = 110_000 # (m)
        return T9 + LK9 * (Z - Z9)
    elseif Z <= 1_000_000
        T10 = 360 # (K)
        Z10 = 120_000 # (m)
        Tinf = 1000 # (K)
        λ = 0.01875 / 1000 # (1 / m)
        ξ = (Z - Z10) * (r0 + Z10) / (r0 + Z)
        return Tinf - (Tinf - T10) * exp(-λ * ξ)
    end
end

# based on David's code, which was based on Regan
# M and log(P) are interpolated quadratically
const Ztableupper = [86000.0, 87000.0, 88000.0, 89000.0, 90000.0, 91000.0,
    93000.0, 95000.0, 97000.0, 99000.0, 101000.0, 103000.0, 105000.0, 107000.0,
    109000.0, 110000.0, 111000.0, 112000.0, 113000.0, 114000.0, 115000.0,
    116000.0, 117000.0, 118000.0, 119000.0, 120000.0, 125000.0, 130000.0,
    135000.0, 140000.0, 145000.0, 150000.0, 160000.0, 170000.0, 180000.0,
    190000.0, 200000.0, 210000.0, 220000.0, 230000.0, 240000.0, 250000.0,
    260000.0, 270000.0, 280000.0, 290000.0, 300000.0, 310000.0, 320000.0,
    330000.0, 340000.0, 350000.0, 360000.0, 370000.0, 380000.0, 390000.0,
    400000.0, 410000.0, 420000.0, 430000.0, 440000.0, 450000.0, 460000.0,
    470000.0, 480000.0, 490000.0, 500000.0, 525000.0, 550000.0, 575000.0,
    600000.0, 625000.0, 650000.0, 675000.0, 700000.0, 725000.0, 750000.0,
    775000.0, 800000.0, 825000.0, 850000.0, 875000.0, 900000.0, 925000.0,
    950000.0, 975000.0, 1000000.0] # (m)
const Ptableupper = [3.7338E-1, 3.1259E-1, 2.6173E-1, 2.1919E-1, 1.8359E-1,
    1.5381E-1, 1.0801E-1, 7.5966E-2, 5.3571E-2, 3.7948E-2, 2.7192E-2,
    1.9742E-2, 1.4477E-2, 1.0751E-2, 8.1142E-3, 7.1042E-3, 6.2614E-3,
    5.5547E-3, 4.9570E-3, 4.4473E-3, 4.0096E-3, 3.6312E-3, 3.3022E-3,
    3.0144E-3, 2.7615E-3, 2.5382E-3, 1.7354E-3, 1.2505E-3, 9.3568E-4,
    7.2028E-4, 5.6691E-4, 4.5422E-4, 3.0395E-4, 2.1210E-4, 1.5271E-4,
    1.1266E-4, 8.4736E-5, 6.4756E-5, 5.0149E-5, 3.9276E-5, 3.1059E-5,
    2.4767E-5, 1.9894E-5, 1.6083E-5, 1.3076E-5, 1.0683E-5, 8.7704E-6,
    7.2285E-6, 5.9796E-6, 4.9630E-6, 4.1320E-6, 3.4498E-6, 2.8878E-6,
    2.4234E-6, 2.0384E-6, 1.7184E-6, 1.4518E-6, 1.2291E-6, 1.0427E-6,
    8.8645E-7, 7.5517E-7, 6.4468E-7, 5.5155E-7, 4.7292E-7, 4.0642E-7,
    3.5011E-7, 3.0236E-7, 2.1200E-7, 1.5137E-7, 1.1028E-7, 8.2130E-8,
    6.2601E-8, 4.8865E-8, 3.9048E-8, 3.1908E-8, 2.6611E-8, 2.2599E-8,
    1.9493E-8, 1.7036E-8, 1.5051E-8, 1.3415E-8, 1.2043E-8, 1.0873E-8,
    9.8635E-9, 8.9816E-9, 8.2043E-9, 7.5138E-9] # (Pa)
const logPtableupper = log.(Ptableupper)
const Mtableupper = [28.95, 28.95, 28.94, 28.93, 28.91, 28.89, 28.82, 28.73,
    28.62, 28.48, 28.30, 28.10, 27.88, 27.64, 27.39, 27.27, 27.14, 27.02,
    26.90, 26.79, 26.68, 26.58, 26.48, 26.38, 26.29, 26.20, 25.80, 25.44,
    25.09, 24.75, 24.42, 24.10, 23.49, 22.90, 22.34, 21.81, 21.30, 20.83,
    20.37, 19.95, 19.56, 19.19, 18.85, 18.53, 18.24, 17.97, 17.73, 17.50,
    17.29, 17.09, 16.91, 16.74, 16.57, 16.42, 16.27, 16.13, 15.98, 15.84,
    15.70, 15.55, 15.40, 15.25, 15.08, 14.91, 14.73, 14.54, 14.33, 13.76,
    13.09, 12.34, 11.51, 10.62, 9.72, 8.83, 8.00, 7.24, 6.58, 6.01, 5.54, 5.16,
    4.85, 4.60, 4.40, 4.25, 4.12, 4.02, 3.94] # (kg / kmol)

function interpolation_index(Z)
    # Find the index for the lower side of the altitude interval
    i = 1
    while i < length(Ztableupper) && Z > Ztableupper[i+1]
        i += 1
    end

    # We are going to reference all elements from i - 1 to i + 1, so we need to
    # adjust the index away from the boundaries
    i == 1 && (i = 2)

    i
end

function interpolation_scale_factors(i, Z)
    Z0 = Ztableupper[i - 1]
    Z1 = Ztableupper[i]
    Z2 = Ztableupper[i + 1]

    scale0 = (Z - Z1) * (Z - Z2) / ((Z0 - Z1) * (Z0 - Z2))
    scale1 = (Z - Z0) * (Z - Z2) / ((Z1 - Z0) * (Z1 - Z2))
    scale2 = (Z - Z0) * (Z - Z1) / ((Z2 - Z0) * (Z2 - Z1))

    scale0, scale1, scale2
end

function pressure_upper(Z)
    i = interpolation_index(Z)
    scale0, scale1, scale2 = interpolation_scale_factors(i, Z)
    logP0 = logPtableupper[i - 1]
    logP1 = logPtableupper[i]
    logP2 = logPtableupper[i + 1]
    logP = logP0 * scale0 + logP1 * scale1 + logP2 * scale2
    exp(logP)
end

function mean_molecular_weight_upper(Z)
    i = interpolation_index(Z)
    scale0, scale1, scale2 = interpolation_scale_factors(i, Z)
    M0 = Mtableupper[i - 1]
    M1 = Mtableupper[i]
    M2 = Mtableupper[i + 1]
    M0 * scale0 + M1 * scale1 + M2 * scale2
end

const speed_of_sound_86km = let
    Z = 86_000.0 # (m)
    H = geopotential_altitude(Z)
    M = mean_molecular_weight_lower(Z)
    T = temperature_lower(H, M)
    speed_of_sound_lower(T, M)
end

immutable State
    mean_molecular_weight::Float64
    temperature::Float64
    pressure::Float64
    speed_of_sound::Float64
end

mean_molecular_weight(s::State) = s.mean_molecular_weight
temperature(s::State) = s.temperature
pressure(s::State) = s.pressure
density(s::State) = pressure(s) * mean_molecular_weight(s) / (Rstar * temperature(s))
speed_of_sound(s::State) = s.speed_of_sound

function check_altitude(Z)
    Z < -5000 && error("altitude below lower bound of -5000 m")
    Z > 1_000_000 && error("altitude above upper bound of 1000000 m")
end

function atmosphere(Z)
    check_altitude(Z)
    if Z < 86_000
        H = geopotential_altitude(Z)
        M = mean_molecular_weight_lower(Z)
        T = temperature_lower(H, M)
        P = pressure_lower(H)
        c = speed_of_sound_lower(T, M)
    else
        T = temperature_upper(Z)
        P = pressure_upper(Z)
        M = mean_molecular_weight_upper(Z)
        c = speed_of_sound_86km
    end
    State(M, T, P, c)
end

end
