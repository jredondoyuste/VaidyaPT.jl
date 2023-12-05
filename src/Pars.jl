module Pars

export read_pars, change_par

include("Background.jl")

import .Background: mass!

"""
    parseable(value): Check if a string can be parsed as an integer.
"""
function parseable(value)
    try
        parse(Int, value)
        return true
    catch
        return false
    end
end

"""
    read_pars(filename): Read the parameters from a file.
"""
function read_pars(filename::String = "pars.txt")
    f = open(filename, "r")
    params = Dict{String,Any}()
    for line in eachline(f)
        if occursin("#", line)
            continue
        end
        parts = split(line, "=")
        name = strip(parts[1])
        value = strip(parts[2])
        if occursin(".", value)
            value = parse(Float64, value)
        elseif parseable(value)
            value = parse(Int, value)
        else
            value = strip(value, ['"'])
        end
        params[name] = value
    end

    params["Δu"] = params["Delta_u"]
    params["Δv"] = params["Delta_v"]
    params["τ"] = params["tau"]
    params["σ_IC"] = params["sigma_IC"]
    params["Matter_σ"] = params["Matter_sigma"]


    close(f)

    params["u_grid"] = collect(range(params["u_0"], params["u_end"], step = params["Δu"]))
    params["v_grid"] = collect(range(params["v_0"], params["v_end"], step = params["Δv"]))
    params["N_u"] = length(params["u_grid"])
    params["N_v"] = length(params["v_grid"])
    params["mass"] =
        mass!(params["v_grid"], params["m_1"], params["m_2"], params["v_1"], params["τ"])

    return params
end

"""
    change_par(p, which, value): Change a parameter in the dictionary `p`.
"""
function change_par(p::Dict{String,Any}, which::String, value::Float64)::Dict{String,Any}
    p[which] = value
    if which in ["u_0", "u_end", "Δu"]
        p["u_grid"] = collect(range(p["u_0"], p["u_end"], step = p["Δu"]))
        p["N_u"] = length(p["u_grid"])
    elseif which in ["v_0", "v_end", "Δv", "m_1", "m_2", "v_1", "τ"]
        p["v_grid"] = collect(range(p["v_0"], p["v_end"], step = p["Δv"]))
        p["N_v"] = length(p["v_grid"])
        p["mass"] = mass!(p["v_grid"], p["m_1"], p["m_2"], p["v_1"], p["τ"])
    end
    if which in ["Delta_u", "Delta_v", "tau", "sigma_IC", "Matter_sigma"]
        p["Δu"] = p["Delta_u"]
        p["Δv"] = p["Delta_v"]
        p["τ"] = p["tau"]
        p["σ_IC"] = p["sigma_IC"]
        p["Matter_σ"] = p["Matter_sigma"]
    end
    return p
end

end