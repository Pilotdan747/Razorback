# Script to run throuugh tutrials for Metaheuristics.jl
# Author: Daniel Owen
# Created: June 28, 2024
# Edited: June 28, 2024
# Version: 0.1

using Metaheuristics

function main()
    f(x) = 10length(x) + sum( x.^2 - 10cos.(2π*x) )

    bounds = BoxConstrainedSpace(lb = -5ones(10), ub = 5ones(10))

    information = Information(f_optimum = 0.0)

    options = Options(f_calls_limit = 9000*10, f_tol = 1e-5)

    algorithm = ECA(information = information, options = options)

    result = optimize(f, bounds, algorithm)

    fx = minimum(result)
    x = minimizer(result)

    x, fx
end

main()
