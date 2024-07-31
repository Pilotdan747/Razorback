# Initial test at propataging planet orbits. Simple circles fornow
# Author: Daniel Owen
# Created: May 14, 2024
# Edited: July 24, 2024

include("kepler.jl")

function test_planet()
    R1, V1 = get_state_circ("Earth", 0)
    R2, V2 = get_state_circ("Earth", 365.25/4*86400)
    R3, V3 = get_state_circ("Earth", 365.25/2*86400)

    R4, V4 = get_state_kep("Earth", 0)
    R5, V5 = get_state_kep("Earth", 365.25)

    # println(R1)
    # println(R2)
    # println(R3)
    println(R4)
    println(V4)
    println(R5)
    println(V5)
end

ids = Dict(1 => "Mercury",
           2 => "Venus",
           3 => "Earth",
           4 => "Mars",
           5 => "Jupiter",
           6 => "Saturn",
           7 => "Uranus",
           8 => "Nepture")

sma = Dict("Mercury" => 0.387098309, 
           "Venus" => 0.72332982, 
           "Earth" => 1.0000010178, 
           "Mars" => 1.52367934,
           "Jupiter" => 5.202603191,
           "Saturn" => 9.554909595,
           "Uranus" => 19.218446061,
           "Nepture" => 30.11038687)

ecc = Dict("Mercury" => 0.205631752, 
           "Venus" => 0.006771882, 
           "Earth" => 0.016708617, 
           "Mars" => 0.093400620,
           "Jupiter" => 0.048494851,
           "Saturn" => 0.055508622,
           "Uranus" => 0.046295898,
           "Nepture" => 0.008988095)
        
inc = Dict("Mercury" => 7.00798625, 
           "Venus" => 3.39446619, 
           "Earth" => 0.000, 
           "Mars" => 1.84972648,
           "Jupiter" => 1.30326966,
           "Saturn" => 2.48887810,
           "Uranus" => 0.77319617,
           "Nepture" => 1.76995221)

RAAN = Dict("Mercury" => 48.33089304, 
            "Venus" => 76.67992019, 
            "Earth" => 0.000, 
            "Mars" => 49.55809321,
            "Jupiter" => 100.46444064,
            "Saturn" => 113.66552370,
            "Uranus" => 74.00594723,
            "Nepture" => 131.78405702)

long_peri = Dict("Mercury" => 77.45611904, 
                 "Venus" => 131.56370724, 
                 "Earth" => 102.93734851, 
                 "Mars" => 336.06023398,
                 "Jupiter" => 14.33130924,
                 "Saturn" => 93.05678728,
                 "Uranus" => 173.00515922,
                 "Nepture" => 48.12369050)

long_true = Dict("Mercury" => 252.25090551, 
                 "Venus" => 181.97980084, 
                 "Earth" => 100.46644851, 
                 "Mars" => 355.43327463,
                 "Jupiter" => 34.35148392,
                 "Saturn" => 50.07747138,
                 "Uranus" => 314.05500511,
                 "Nepture" => 304.34866548)

plan_μ = Dict("Mercury" => 2.2032e4, 
              "Venus" => 3.257e5, 
              "Earth" => 3.986004415e5, 
              "Mars" => 4.305e4, 
              "Jupiter" => 1.268e8, 
              "Saturn" => 3.794e7, 
              "Uranus" => 5.794e6, 
              "Neptune" => 6.809e6)

R_eq = Dict("Mercury" => 2439.0, 
            "Venus" => 6052.0, 
            "Earth" => 6378.1363, 
            "Mars" => 3397.2, 
            "Jupiter" => 71492.0, 
            "Saturn" => 60268.0, 
            "Uranus" => 25559.0, 
            "Neptune" => 24764.0)

μ_sun = 1.32712440e11;
AU = 149597870.7

function get_state_kep(planet, time)
    a = sma[planet]*AU
    e = ecc[planet]
    i = deg2rad(inc[planet])
    Ω = deg2rad(RAAN[planet])
    ω_true = deg2rad(long_peri[planet])
    ω = ω_true - Ω
    λ_true = deg2rad(long_true[planet])
    θ = λ_true - Ω - ω
    u = ω + θ

    special_case = "none"
    if e == 0
        if i == 0
            special_case = "circ_eq"
        else
            special_case = "circ_inc"
        end
    else
        if i == 0
            special_case = "ellip_eq"
        end
    end
    
    coe = Dict("a" => a, 
               "e" => e, 
               "i" => i, 
               "Ω" => Ω, 
               "ω" => ω, 
               "θ" => θ,
               "ω true" => ω_true,
               "λ true" => λ_true,
               "u" => u, 
               "special_case" => special_case)

    R_j2000, V_j2000 = elements_to_state(coe, μ_sun)
    return propogate(R_j2000, V_j2000, μ_sun, time)
end

function get_state_circ(planet, time)
    r = sma[planet]*AU
    v = sqrt(μ_sun/r)

    T = 2*π/sqrt(μ_sun)*r^(3/2)

    θ = 2*π*time/T

    R = [r*cos(θ), r*sin(θ), 0]
    V = [-v*sin(θ), v*cos(θ), 0]

    return (R, V)
end

# test_planet()