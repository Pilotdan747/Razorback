# Initial test at seting up a simple mga itinerary - EVEEJU
# Author: Daniel Owen
# Created: May 14, 2024
# Edited: May 14, 2024
# Version: 0.1

include("lambert.jl")
include("planets.jl")

dt1 = 30*86400
dt2 = 30*86400
dt3 = 900*86400
dt4 = 1400*86400

t_earth = 0
t_venus = t_earth + dt1
t_earth_2 = t_venus + dt2
t_jupiter = t_earth_2 + dt3
t_uranus = t_jupiter + dt4

R_Earth, V_Earth = get_state("Earth", t_earth)
R_Venus, V_Venus = get_state("Venus", t_venus)
R_Earth_2, V_Earth_2 = get_state("Earth", t_earth_2)
R_Jupiter, V_Jupiter = get_state("Jupiter", t_jupiter)
R_Uranus, V_Uranus = get_state("Uranus", t_uranus)

V1, V2 = lambert_battin(R_Earth, R_Venus, dt1, μ_sun, 0)
V3, V4 = lambert_battin(R_Venus, R_Earth, dt2, μ_sun, 0)
V5, V6 = lambert_battin(R_Earth, R_Jupiter, dt3, μ_sun, 0)
V7, V8 = lambert_battin(R_Jupiter, R_Uranus, dt4, μ_sun, 0)

function cost()
    V_∞_Earth_out = norm(V1 - V_Earth)
    V_∞_Venus_in = norm(V2 - V_Venus)
    V_∞_Venus_out = norm(V3 - V_Venus)
    V_∞_Earth_in_2 = norm(V4 - V_Earth_2)
    V_∞_Earth_out_2 = norm(V5 - V_Earth_2)
    V_∞_Jupiter_in = norm(V6 - V_Jupiter)
    V_∞_Jupiter_out = norm(V7 - V_Jupiter)
    V_∞_Uranus_in = norm(V8 - V_Uranus)

    j_balistic = abs(V_∞_Venus_out - V_∞_Venus_in) + abs(V_∞_Earth_out_2 - V_∞_Earth_in_2) + abs(V_∞_Jupiter_out - V_∞_Jupiter_in)
    j_depart = V_∞_Earth_out
    j_arrive = V_∞_Uranus_in

    j = j_balistic + j_depart + j_arrive
    println(j)
end
