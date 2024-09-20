# Early implementation of PSO optimizer on MGA problem
# Author: Daniel Owen
# Created: July 5, 2024
# Edited: July 30, 2024

using Pkg
Pkg.activate(pwd())

using Metaheuristics

include("../razorback/mga.jl")
include("plot_traj.jl")

function optum_wrapper(x)
    init_t = x[1]
    dt_vals = x[2:5]
    
    planets = ["Earth"]
    for id in x[6:end]
        planet = ids[Int(round(id))]
        push!(planets, planet)
    end
    push!(planets, "Jupiter")
    
    data = generate_mga_data(init_t, dt_vals, planets)

    return cost(data, false)
end

function optum_wrapper_parallel(x)
    f = zeros(size(x, 1))
    Threads.@threads for i in i:size(x, 1)
        init_t = x[i, 1]
        dt_vals = x[i, 2:5]
        
        planets = ["Earth"]
        for id in x[i, 6:end]
            planet = ids[Int(round(id))]
            push!(planets, planet)
        end
        push!(planets, "Jupiter")
        
        data = generate_mga_data(init_t, dt_vals, planets)

        f[i] = cost(data, false)
    end

    return f
end

function min_data(x)
    init_t = x[1]
    dt_vals = x[2:5]
    
    planets = ["Earth"]
    for x in x[6:end]
        planet = ids[Int(round(x))]
        push!(planets, planet)
    end
    push!(planets, "Jupiter")

    return generate_mga_data(init_t, dt_vals, planets)
end

function main()
    low_bounds = [25*365.25*86400; 30*86400*ones(4); 2*ones(3)]
    up_bounds = [35*365.25*86400; 10000*86400*ones(4); 5*ones(3)]

    bounds = boxconstraints(lb = low_bounds, ub = up_bounds)
    # println(bounds)

    options = Options(iterations = 50, f_calls_limit = 1e7, debug = true)
    alg = PSO(N = 2500, options = options)
    # alg = SA(N = 2500, options = options)
    # alg = DE(N = 1000, options = options)
    # alg = ABC(N = 2500, options = options)
    # alg = GA(N = 2500) # Broken
    # alg = ECA(N = 100, options = options)


    # # f(x) = x^2 + sin(x)
    # if Threads.nthreads() > 1
    #     result = optimize(optum_wrapper_parallel, bounds, alg)
    # else
    #     result = optimize(optum_wrapper, bounds, alg)
    # end

    result = optimize(optum_wrapper, bounds, alg)

    println(result)

    x = minimizer(result)
    flyby1 = Int(round(x[6]))
    flyby2 = Int(round(x[7]))
    flyby3 = Int(round(x[8]))
    println("Flybys: $flyby1, $flyby2, $flyby3")

    data = min_data(x)
    cost(data, true)

    save_name = ""
    for planet in data.itinerary
        if planet == data.itinerary[end]
            save_name = save_name + "$planet"
        else
            save_name = save_name + "$planet _"
        end
    end

    save_mga_data(data, save_name)

    fig = plot()

    for planet in unique(data.itinerary)
        plot_planet!(fig, planet)
    end



    plot_solution!(fig, data)

    relayout!(fig, scene_aspectmode="data")
    savefig(fig, "mga_traj.png")
    
    html_file = open("mga_traj.html", "w")
    PlotlyBase.to_html(html_file, fig)
    close(html_file)
end

main()