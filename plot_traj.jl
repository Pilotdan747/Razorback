# Initial test at plotting a solution to the MGA problem
# Author: Daniel Owen
# Created: May 23, 2024
# Edited: May 23, 2024
# Version: 0.1

# Prop planets and plot
# Load inital state of each arc
# Kepler to prop arc
# Plot arcs

include("planets.jl")

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

solution = []
arc["dT"] = 0
arc["R"] = []
arc["V"] = []

solution[1] = arc

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
