# Initial test at seting up a simple mga itinerary - EVEEJU
# Author: Daniel Owen
# Created: May 14, 2024
# Edited: July 5, 2024

include("../razorback/lambert.jl")
include("../razorback/planets.jl")

struct mga_data
    position_vectors
    planet_velocities
    traj_velocities
end

function init_state()
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

    pos_vecs = [R_Earth, R_Venus, R_Earth_2, R_Jupiter, R_Uranus]
    plan_vels = [V_Earth, V_Venus, V_Earth_2, V_Jupiter, V_Uranus]
    traj_vels = [V1, V2, V3, V4, V5, V6, V7, V8]

    return mga_data(pos_vecs, plan_vels, traj_vels)
end

function init_state_generic(dt_vals, planet)
    times = [0]
    for (i, dt) in enumerate(dt_vals)
        t = times[i] + dt
        push!(times, t)
    end

    pos_vecs = []
    plan_vels = []
    for (i, time) in enumerate(times)
        R, V = get_state(planet[i], time)
        push!(pos_vecs, R)
        push!(plan_vels, V)
    end

    traj_vels = []
    for (i, dt) in enumerate(dt_vals)
        V1, V2 = lambert_battin(pos_vecs[i], pos_vecs[i + 1], dt, μ_sun, 0)
        push!(traj_vels, V1)
        push!(traj_vels, V2)
    end
    
    return mga_data(pos_vecs, plan_vels, traj_vels)
end

function cost(data)
    V1 = data.traj_velocities[1]
    V2 = data.traj_velocities[2]
    V3 = data.traj_velocities[3]
    V4 = data.traj_velocities[4]
    V5 = data.traj_velocities[5]
    V6 = data.traj_velocities[6]
    V7 = data.traj_velocities[7]
    V8 = data.traj_velocities[8]

    V_Earth = data.planet_velocities[1]
    V_Venus = data.planet_velocities[2]
    V_Earth_2 = data.planet_velocities[3]
    V_Jupiter = data.planet_velocities[4]
    V_Uranus = data.planet_velocities[5]

    V_∞_Earth_out = norm(V1 - V_Earth) # 1
    V_∞_Venus_in = norm(V2 - V_Venus) # 2
    V_∞_Venus_out = norm(V3 - V_Venus) # 2
    V_∞_Earth_in_2 = norm(V4 - V_Earth_2) # 3
    V_∞_Earth_out_2 = norm(V5 - V_Earth_2) # 3
    V_∞_Jupiter_in = norm(V6 - V_Jupiter) # 4
    V_∞_Jupiter_out = norm(V7 - V_Jupiter) # 4
    V_∞_Uranus_in = norm(V8 - V_Uranus) # 5

    j_balistic = abs(V_∞_Venus_out - V_∞_Venus_in) + abs(V_∞_Earth_out_2 - V_∞_Earth_in_2) + abs(V_∞_Jupiter_out - V_∞_Jupiter_in)
    j_depart = V_∞_Earth_out
    j_arrive = V_∞_Uranus_in

    j = j_balistic + j_depart + j_arrive
    println("Cost: $j")
end

function cost_generic(data)
    V∞_arr = []
    
    for i in 2:length(data.traj_velocities) + 1
        traj_idx = i - 1
        planet_idx = convert(Int32, (floor((i + 1)/2)))

        V∞ = norm(data.traj_velocities[traj_idx] - data.planet_velocities[planet_idx])
        push!(V∞_arr, V∞)
    end

    j_balistic = 0
    for i in 2:2:length(data.traj_velocities) -1
        j_balistic += abs(V∞_arr[i+1] - V∞_arr[i])
    end

    j_depart = V∞_arr[1]
    j_arrive = V∞_arr[end]

    j = j_balistic + j_depart + j_arrive
    println("Cost: $j")
end

dt_vals = [30*86400, 30*86400, 900*86400, 1400*86400]
planets = ["Earth", "Venus", "Earth", "Jupiter", "Uranus"]

data = init_state()
cost(data)
data_gen = init_state_generic(dt_vals, planets)
cost_generic(data)
