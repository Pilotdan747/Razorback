# Initial test at propataging planet orbits. Simple circles fornow
# Author: Daniel Owen
# Created: May 14, 2024
# Edited: July 5, 2024
# Version: 0.2

function test_planet()
    R1, V1 = get_state("Earth", 0)
    R2, V2 = get_state("Earth", 365.25/4*86400)
    R3, V3 = get_state("Earth", 365.25/2*86400)

    println(R1)
    println(R2)
    println(R3)
end

ids = Dict(1 => "Mercury",
           2 => "Venus",
           3 => "Earth",
           4 => "Mars",
           5 => "Jupiter",
           6 => "Saturn",
           7 => "Uranus",
           8 => "Nepture")

sma = Dict("Mercury" => 0.38709927, 
           "Venus" => 0.72333566, 
           "Earth" => 1.00000261, 
           "Mars" => 1.52371034,
           "Jupiter" => 5.20288700,
           "Saturn" => 9.53667594,
           "Uranus" => 19.18916464,
           "Nepture" => 30.06992276)

μ_sun = 1.32712440e11;
AU = 149597870.7

function get_state(planet, time)
    r = sma[planet]*AU
    v = sqrt(μ_sun/r)

    T = 2*π/sqrt(μ_sun)*r^(3/2)

    θ = 2*π*time/T

    R = [r*cos(θ), r*sin(θ), 0]
    V = [-v*sin(θ), v*cos(θ), 0]

    return (R, V)
end