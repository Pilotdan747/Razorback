# Initial test at plotting a solution to the MGA problem
# Author: Daniel Owen
# Created: May 23, 2024
# Edited: May 23, 2024
# Version: 0.1

# Prop planets and plot
# Load inital state of each arc
# Kepler to prop arc
# Plot arcs

using PlotlyJS

include("../razorback/planets.jl")
include("../razorback/kepler.jl")

function plot_planet!(plot, planet)
    a = sma[planet]
    T = 2*π/sqrt(μ_sun)*a^(3/2)

    one_rev = range(0, T, step=86400)

    states = []
    for time in one_rev
        r, v = get_state(planet, time)
        push!(states, r)
    end

    # add each state to a 3d plot
    # I think plotly needs it as a trace object?

end

# solution = []
# arc["dT"] = 0
# arc["R"] = []
# arc["V"] = []

# solution[1] = arc

function plot_solution!(plot, solution)
    for arc in solution
        r = solution["R"][i]
        v = solution["V"][i]
        dT = solution["dT"][i]

        time_range = range(0, 86400, dT)
        for time in time_range
            # kep solve with x0 = r, v for dT seconds
            # add R's to trace
        end
         
        # add trace to plot
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

simple_plot()