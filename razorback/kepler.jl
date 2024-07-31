# Initial test at basic kepler's problem functions (prop and state <=> element)
# Author: Daniel Owen
# Created: May 23, 2024
# Edited: July 24, 2024
# Version: 0.2

using LinearAlgebra

function state_to_elements(R, V, μ)
    H = cross(R, V)
    h = norm(H)

    N = cross([0, 0, 1], H)
    n = norm(N)

    r = norm(R)
    v = norm(V)

    ϵ = v^2/2 - μ/r
    E = ((v^2 - μ/r)*R - dot(R, V)*V)/μ
    e = norm(E)

    a = -μ/2/ϵ

    i = acos(H[3]/h)
    
    Ω = acos(N[1]/n)
    if N[2] < 0
        Ω = 2*π - Ω
    end

    ω = acos(dot(N, E)/n/e)
    if E[3] < 0
        ω = 2*π - ω
    end

    try
        θ = acos(dot(E, R)/e/r)
    catch DomainError
        θ = acos(1.0)
    end

    if dot(R, V) < 0
        θ = 2*π - θ
    end

    special_case = "none"
    λ_true = 0
    u = 0
    ω_true = 0
    if e == 0
        if i == 0
            special_case = "circ_eq"
            λ_true = acos(R[1]/r)
            if R[2] < 0
                λ_true = 2*π - λ_true
            end
        else
            special_case = "circ_inc"
            u = acos(dot(N, R)/n/r)
            if R[3] < 0
                u = 2*π - u
            end
        end
    else
        if i == 0
            special_case = "ellip_eq"
            ω_true = acos(E[1]/e)
            if E[2] < 0
                ω_true = 2*π - ω_true
            end
        end
    end

    OE = Dict("a" => a, 
              "e" => e, 
              "i" => i, 
              "Ω" => Ω, 
              "ω" => ω, 
              "θ" => θ,
              "ω true" => ω_true,
              "λ true" => λ_true,
              "u" => u, 
              "special_case" => special_case)
    return OE
end

function elements_to_state(OE, μ)
    θ = OE["θ"]
    Ω = OE["Ω"]
    ω = OE["ω"]

    special_case = OE["special_case"]
    if special_case != "none"
        if special_case == "circ_eq"
            Ω = 0
            ω = 0
            θ = OE["λ true"]
        elseif special_case == "circ_inc"
            ω = 0
            θ = OE["u"]
        else
            Ω = 0
            ω = OE["ω true"]
        end
    end

    h = sqrt(μ*OE["a"]*(1 - OE["e"]^2))
    r = h^2/μ/(1 + OE["e"]*cos(θ))
    
    ϵ = -μ/2/OE["a"]
    v = sqrt(2*(ϵ + μ/r))

    # println(ϵ)
    # println(v)

    R_perifocal = [r*cos(θ), r*sin(θ), 0]
    V_perifocal = [v*-sin(θ), v*cos(θ), 0]
    
    # rotate w/ matrix
    R1 = [cos(Ω) sin(Ω) 0; -sin(Ω) cos(Ω) 0; 0 0 1]'
    R2 = [1 0 0; 0 cos(OE["i"]) sin(OE["i"]); 0 -sin(OE["i"]) cos(OE["i"])]'
    R3 = [cos(ω) sin(ω) 0; -sin(ω) cos(ω) 0; 0 0 1]'

    R = R1*R2*R3*R_perifocal
    V = R1*R2*R3*V_perifocal

    return (R, V)
end

function kep_E(M, e)
    if M > π || -π < M < 0 
        E_old = M - e
    else
        E_old = M + e
    end

    for i in 1:100
        E = E_old + (M - E_old + e*sin(E_old))/(1 - e*cos(E_old))
        
        if abs(E - E_old) < 1e-9
            break
        end
        
        E_old = E
    end

    return E_old
end

function kep_H(M, e)
end

function propogate(R, V, μ, dT)
    OE = state_to_elements(R, V, μ)
    e = OE["e"]

    if OE["a"] < 0
        n = 2*sqrt(μ/(-1*OE["a"]^3))
        Hi = acosh((e + cos(OE["θ"]))/(1 + e*cos(OE["θ"])))
        Mi = e*sinh(Hi) - Hi
    else
        n = 2*sqrt(μ/OE["a"]^3)
        Ei = 2*atan(sqrt((1 - e)/(1 + e))*tan(OE["θ"]/2))
        Mi = Ei - e*sin(Ei)
    end

    Mf = Mi + dT*n
    Mf = mod(Mf, 2*π)

    if OE["a"] < 0
        Hf = kep_H(Mf, e)

        OE["θ"] = 1
    else
        Ef = kep_E(Mf, e)

        OE["θ"] = 2*atan(sqrt((1 + e)/(1 - e))*tan(Ef/2))
    end

    return elements_to_state(OE, μ)
end