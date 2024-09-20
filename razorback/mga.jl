# Library functions for setting up MGA optimization problem
# Author: Daniel Owen
# Created: July 5, 2024
# Edited: July 24, 2024

using JSON

include("lambert.jl")
include("planets.jl")

function turning_angle(V_in, V_out, planet)
    Re = R_eq[planet]
    μ = plan_μ[planet]

    turn = 0
    try
        turn = acos(dot(V_in, V_out)/norm(V_in)/norm(V_out))
    catch DomainError
        turn = acos(1.0)
    end

    max_turn = 2*asin(1/(1 + Re*norm(V_in)^2/μ))

    angle_diff = turn - max_turn
    return max(angle_diff, 0.0)
end

function cost(data, print)
    V∞_arr = []
    
    for i in 2:length(data.traj_velocities) + 1
        traj_idx = i - 1
        planet_idx = convert(Int32, (floor((i + 1)/2)))

        V∞ = norm(data.traj_velocities[traj_idx] - data.planet_velocities[planet_idx])
        push!(V∞_arr, V∞)
    end

    j_balistic = 0
    j_turn = 0
    for i in 2:2:length(data.traj_velocities) - 1
        planet_idx = convert(Int32, (floor((i + 1)/2)))
        planet = data.itinerary[planet_idx]

        j_balistic += abs(V∞_arr[i+1] - V∞_arr[i])
        j_turn += turning_angle(V∞_arr[i+1], V∞_arr[i], planet)
    end

    j_depart = V∞_arr[1]
    j_arrive = V∞_arr[end]

    j = 100*j_balistic + 100*j_turn + j_depart + j_arrive

    if print
        println("J: $j   Ballistic: $j_balistic   Turn: $j_turn     Departure: $j_depart   Arrival: $j_arrive")
    end
    return j
end

struct mga_data
    position_vectors
    planet_velocities
    traj_velocities
    itinerary
    dt_vals
end

function save_mga_data(data, save_name)
    j_str = JSON.json(data)

    # TODO how to specify file?
    save_str = "output/$save_name"
    f = open("mga_problem.json", "w")
    write(f, j_str)
    close(f)
end

function load_mga_data(file_str)
    json_dict = JSON.parse(readline(file_str))

    data = mga_data(json_dict["position_vectors"],
                    json_dict["planet_velocities"],
                    json_dict["traj_velocities"],
                    json_dict["itinerary"],
                    json_dict["dt_vals"])

    return data
end

function generate_mga_data(init_time, dt_vals, planet)
    times = [init_time]
    for (i, dt) in enumerate(dt_vals)
        t = times[i] + dt
        push!(times, t)
    end

    pos_vecs = []
    plan_vels = []
    for (i, time) in enumerate(times)
        # R, V = get_state_circ(planet[i], time)
        R, V = get_state_kep(planet[i], time)
        push!(pos_vecs, R)
        push!(plan_vels, V)
    end

    traj_vels = []
    for (i, dt) in enumerate(dt_vals)
        V1, V2 = lambert_battin(pos_vecs[i], pos_vecs[i + 1], dt, μ_sun, 0)
        push!(traj_vels, V1)
        push!(traj_vels, V2)
    end
    
    return mga_data(pos_vecs, plan_vels, traj_vels, planet, dt_vals)
end