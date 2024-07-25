# Early implementation of PSO optimizer on MGA problem
# Author: Daniel Owen
# Created: July 5, 2024
# Edited: July 23, 2024
# Version: 0.1

using Metaheuristics

include("../razorback/mga.jl")

function optum_wrapper(x)
    init_t = x[1]
    dt_vals = x[2:5]
    
    planets = ["Earth"]
    for id in x[6:end]
        planet = ids[Int(round(id))]
        push!(planets, planet)
    end
    push!(planets, "Uranus")
    
    data = generate_mga_data(init_t, dt_vals, planets)

    return cost(data, false)
end

function min_data(x)
    init_t = x[1]
    dt_vals = x[2:5]
    
    planets = ["Earth"]
    for x in x[6:end]
        planet = ids[Int(round(x))]
        push!(planets, planet)
    end
    push!(planets, "Uranus")

    return generate_mga_data(init_t, dt_vals, planets)
end

function main()
    low_bounds = [25*365.25*86400; 30*86400*ones(4); 2*ones(3)]
    up_bounds = [35*365.25*86400; 10000*86400*ones(4); 7*ones(3)]

    bounds = boxconstraints(lb = low_bounds, ub = up_bounds)
    # println(bounds)

    options = Options(iterations = 1000, f_calls_limit = 1e7, debug = true)
    # alg = PSO(N = 100, options = options)
    # alg = SA(x_initial = (up_bounds + low_bounds)/2, N = 1000, options = options)
    # alg = DE(N = 1000, options = options)
    # alg = ABC(N = 2500, options = options)
    # alg = GA(N = 2500) # Broken
    alg = ECA(N = 100, options = options)


    # f(x) = x^2 + sin(x)
    result = optimize(optum_wrapper, bounds, alg)
    println(result)

    x = minimizer(result)
    flyby1 = Int(round(x[6]))
    flyby2 = Int(round(x[7]))
    flyby3 = Int(round(x[8]))
    println("Flybys: $flyby1, $flyby2, $flyby3")

    data = min_data(x)
    cost(data, true)
end

main();