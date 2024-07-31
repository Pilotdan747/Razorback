# Initial test at plotting a solution to the MGA problem
# Author: Daniel Owen
# Created: May 23, 2024
# Edited: July 24, 2024
# Version: 0.1


using PlotlyJS

include("../razorback/planets.jl")
include("../razorback/kepler.jl")

function plot_planet!(fig, planet)
    a = sma[planet]*AU
    T = 2*π/sqrt(μ_sun)*a^(3/2)

    one_rev = range(0, T, step=86400)

    X = []
    Y = []
    Z = []
    for time in one_rev
        r, v = get_state_kep(planet, time)
        push!(X, r[1]/AU)
        push!(Y, r[2]/AU)
        push!(Z, r[3]/AU)
    end

    plan_trace = scatter3d(x=X, y=Y, z=Z, mode="lines", name=planet, aspect_ratio = 1)
    add_trace!(fig, plan_trace)
end

function plot_solution!(fig, solution)
    itinerary = solution.itinerary
    pos_vec = solution.position_vectors
    vel_vec = solution.traj_velocities
    dt_vals = solution.dt_vals

    for i = 1:length(itinerary) - 1
        r = pos_vec[i]
        v = vel_vec[i]
        dT = dt_vals[i]

        planet = itinerary[i]
        next_planet = itinerary[i + 1]


        name = "$planet to $next_planet"
        println(name)

        time_range = range(0, dT, step=86400)

        X = []
        Y = []
        Z = []
        for time in time_range
            r_t, v_t = propogate(r, v, μ_sun, time)
            push!(X, r_t[1]/AU)
            push!(Y, r_t[2]/AU)
            push!(Z, r_t[3]/AU)
        end
         
        trace = scatter3d(x = X, y = Y, z = Z, mode = "lines", name = name)
        add_trace!(fig, trace)
    end
end

function simple_plot()
    N = 1000
    t_vals = range(0, 6*3600, N)
    
    X_arr = []
    Y_arr = []
    Z_arr = []
    X_arr2 = []
    Y_arr2 = []
    Z_arr2 = []
    OE = Dict("a" => 7000, "e" => 0.8, "i" => π/4, "Ω" => π/4, "ω" => π/4, "θ" => 0, "special_case" => "none")
    OE2 = Dict("a" => 7000, "e" => 0.2, "i" => π/2, "Ω" => π/4, "ω" => π/4, "θ" => 0, "special_case" => "none")
    R, V = elements_to_state(OE, 3.986e5)
    R2, V2 = elements_to_state(OE2, 3.986e5)
    for t in t_vals
        R_t, V_t = propogate(R, V, 3.986e5, t)
        R_t2, V_t2 = propogate(R2, V2, 3.986e5, t)
        
        push!(X_arr, R_t[1])
        push!(Y_arr, R_t[2])
        push!(Z_arr, R_t[3])

        push!(X_arr2, R_t2[1])
        push!(Y_arr2, R_t2[2])
        push!(Z_arr2, R_t2[3])
    end

    trace1 = scatter3d(x=X_arr, y=Y_arr, z=Z_arr, mode="lines", name="Orbit 1")
    trace2 = scatter3d(x=X_arr2, y=Y_arr2, z=Z_arr2, mode="lines", name="Orbit 2")

    fig = plot([trace1, trace2])

    # open("test.html", "w") do f
    #     PlotlyBase.to_html(f, fig.plot)
    # end
end

# simple_plot()
fig = plot()
planets = ["Venus", "Earth", "Jupiter"]

for p in planets
    plot_planet!(fig, p)
end

relayout!(fig, scene_aspectmode="data")
fig