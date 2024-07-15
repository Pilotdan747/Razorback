# Library functions for setting up MGA optimization problem
# Author: Daniel Owen
# Created: July 5, 2024
# Edited: July 6, 2024
# Version: 0.1

include("lambert.jl")
include("planets.jl")

function cost(data, print)
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

    if print
        println("J: $j   Ballistic: $j_balistic   Departure: $j_depart   Arrival: $j_arrive")
    end
    return j
end

struct mga_data
    position_vectors
    planet_velocities
    traj_velocities
end

function generate_mga_data(dt_vals, planet)
    times = [0.0]
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