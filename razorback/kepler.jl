# Initial test at basic kepler's problem functions (prop and state <=> element)
# Author: Daniel Owen
# Created: May 23, 2024
# Edited: May 23, 2024
# Version: 0.1

function state_to_elements(R, V, μ)
    # calc h
    # calc v
    # calc a from energy
    # e from h equation
    # look up i, raan, argPeri
end

function elements_to_state(OE, μ)
    h = sqrt(mu*OE["a"]*(1 - OE["e"]^2))
    r = h^2/μ/(1 + OE["e"]*cos(OE["θ"]))
    
    ϵ = -μ/2/OE["a"]
    v = sqrt(2*(ϵ + μ/r))

    R = [r*cos(OE["θ"]), r*sin(OE["θ"]), 0]
    V = [v*-sin(OE["θ"]), v*cos(OE["θ"]), 0]
    
    # rotate w/ matrix
end

function propogate(R, V, μ, dT)
    OE = state_to_elements(R, V, μ)
    T = 2*π/sqrt(μ)*OE["a"]^(3/2)

    # get Ei
    Mi = Ei - OE["e"]*sin(Ei)
    Mf = Mi + dT*2*pi/T

    # get new_theta

    OE["θ"] = new_theta
    
    return elements_to_state(OE, μ)
end