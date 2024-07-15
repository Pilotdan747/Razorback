# Early implementation of PSO optimizer on MGA problem
# Author: Daniel Owen
# Created: July 5, 2024
# Edited: July 5, 2024
# Version: 0.1

using Metaheuristics

include("../razorback/mga.jl")

function optum_wrapper(x)
    dt_vals = x[1:4]
    
    planets = ["Earth"]
    for x in x[5:end]
        planet = ids[Int(round(x))]
        push!(planets, planet)
    end
    
    data = generate_mga_data(dt_vals, planets)

    return cost(data, false)
end

function min_data(x)
    dt_vals = x[1:4]
    
    planets = ["Earth"]
    for x in x[5:end]
        planet = ids[Int(round(x))]
        push!(planets, planet)
    end

    return generate_mga_data(dt_vals, planets)
end

function main()
    low_bounds = [30*86400*ones(4); 2*ones(4)]
    up_bounds = [1500*86400*ones(4); 8*ones(4)]

    bounds = boxconstraints(lb = low_bounds, ub = up_bounds)
    # println(bounds)

    options = Options(iterations = 2500, f_calls_limit = 1000000)
    alg = PSO(N = 250, options = options)

    # f(x) = x^2 + sin(x)
    result = optimize(optum_wrapper, bounds, alg)
    println(result)

    x = minimizer(result)

    data = min_data(x)
    cost(data, true)
end

main()