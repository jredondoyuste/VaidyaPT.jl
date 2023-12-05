module VaidyaPT

export Field, Initialize, Update, Evolve, read_pars, change_par, mismatch, minimize, find_N_cycles, SimData, Inference_Config, inference, read_data, shift, cut, shift_cut, reconstruct, get_vshift, add_noise

include("Pars.jl")
include("PerturbationEquations.jl")
include("Background.jl")
include("Minimization.jl")
include("Nested_Sampling.jl")

import .Pars: read_pars, change_par
import .PerturbationEquations: Source!, Potential!
import .Background: r_sol!, f_sol!
import .Minimization: mismatch, minimize, find_N_cycles
import .Nested_Sampling: SimData, Inference_Config, inference, read_data, shift, cut, shift_cut, reconstruct, get_vshift, add_noise

using StructArrays
using Parameters
using Term
using ProgressMeter

"""
    gaussian_profile!(A, v, v0, σ): Gaussian profile.
"""
function gaussian_profile!(
    A::Float64,
    v::Union{Float64,Array{Float64}},
    v0::Float64,
    σ::Float64,
    )::Union{Float64,Array{Float64}}
    return @. A * exp(-(v - v0)^2 / (2 * σ^2))
end

"""
    Field(parfile): Field structure.
"""
@with_kw mutable struct Field
    parfile::String = "pars.txt"
    p::Dict{String,Any} = read_pars(parfile)
    ψ::Array{Float64,1} = zeros(Float64, p["N_v"])
    u::Float64 = p["u_0"]
    ϕ::Array{Float64,1} = zeros(Float64, p["N_v"])
end

"""
    Initialize(U): Initialize the field.
"""
function Initialize(U::Field)
    U.ψ = gaussian_profile!(U.p["A_IC"], U.p["v_grid"], U.p["v_IC"], U.p["σ_IC"])
    U.ϕ = gaussian_profile!(U.p["Matter_A"], U.p["v_grid"], U.p["Matter_v0"], U.p["Matter_σ"])
    #--------------------------------------
    tprint(
        Panel(
            "\n",
            Panel(
                "A = $(U.p["A_IC"]), v_1 = $(U.p["v_IC"]), sigma = $(U.p["sigma_IC"])",
                title = "Field",
                title_style = "bold",
                style = "blue",
                box = :ROUNDED,
                fit = true,
            ),
            Panel(
                "A = $(U.p["Matter_A"]), v_1 = $(U.p["Matter_v0"]), sigma = $(U.p["Matter_sigma"])",
                title = "Matter",
                title_style = "bold",
                style = "blue",
                box = :ROUNDED,
                fit = true,
            ),
            title = "Parameters",
            title_style = "bold",
            style = "bold red",
            box = :ROUNDED,
            justify = :center,
            fit = true,
        ),
    )
    #--------------------------------------
    return nothing
end

"""
    Update(U): Update the field one time step.
"""
function Update(U::Field)::Field
    U_new = Field(
        p = U.p,
        ψ = zeros(Float64, U.p["N_v"]),
        u = U.u + U.p["Δu"],
        ϕ = U.ϕ,
    )

    U_new.ψ[1] = 0

    r_vals = r_sol!(
        U.u,
        U.p["m_1"],
        U.p["m_2"],
        U.p["v_0"],
        U.p["v_1"],
        U.p["v_end"],
        U.p["tau"],
        U.p["Δv"],
    )
    f_vals = f_sol!(
        U.u,
        U.p["m_1"],
        U.p["m_2"],
        U.p["v_0"],
        U.p["v_1"],
        U.p["v_end"],
        U.p["tau"],
        U.p["Δv"],
        U.p["Δu"],
    )
    src = Source!(r_vals, f_vals, U.ϕ)

    for i = 2:U.p["N_v"]
        V = f_vals[i] * Potential!(U.p["s"], U.p["l"], U.p["mass"][i], r_vals[i])
        U_new.ψ[i] = U.ψ[i] + U_new.ψ[i-1] - U.ψ[i-1] -
            (U.p["Δv"] * U.p["Δu"] / 2) * V * (U.ψ[i] + U_new.ψ[i-1]) -
            (U.p["Δv"] * U.p["Δu"] / 2) * src[i-1]
    end

    return U_new
end

"""
    Evolve(U, directory_save): Evolve the field and save the data in the folder directory_save.
"""
function Evolve(U::Field, directory_save::String)
    #--------------------------------------
    #--------------------------------------
    tprint(
        Panel(
            "Initializing Evolution...",
            title = "Evolution",
            title_style = "bold",
            box = :ROUNDED,
            fit = true,
        ),
    )
    #--------------------------------------
    tprint(
        Panel(
            Panel(
                "m_1 = $(U.p["m_1"]), m_2 = $(U.p["m_2"]), v_1 = $(U.p["v_1"]), tau = $(U.p["tau"])",
                title = "Background Mass Profile",
                style = "blue",
                title_style = "bold blue",
                fit = true,
            ),
            Panel(
                "u_0 = $(U.p["u_0"]), u_end = $(U.p["u_end"]), Δu = $(U.p["Δu"]), \nv_0 = $(U.p["v_0"]), v_end = $(U.p["v_end"]), Δv = $(U.p["Δv"])",
                title = "Grid",
                style = "blue",
                title_style = "bold blue",
                justify = :center,
                fit = true,
            ),
            box = :ROUNDED,
            style = "red",
            fit = true,
        ),
    )
    #--------------------------------------
    #--------------------------------------
    filename_H =
        directory_save *
        "/" *
        "Horizon_l_" *
        string(U.p["l"]) *
        "_s_" *
        string(U.p["s"]) *
        ".dat"
    filename_I =
        directory_save *
        "/" *
        "Scri_l_" *
        string(U.p["l"]) *
        "_s_" *
        string(U.p["s"]) *
        ".dat"
    filename_R = 
        directory_save *
        "/" *
        "R_"  * string(U.p["r_ext"]) * 
        "_l_" * string(U.p["l"]) *
        "_s_" * string(U.p["s"]) *
        ".dat"
    #--------------------------------------
    open(directory_save * "/pars.txt", "w") do f
        for (key, value) in U.p
            if key in ["u_grid", "v_grid", "mass"]
                continue
            end
            write(f, string(key, " = ", value, "\n"))
        end
    end
    #--------------------------------------
    ψ_scri = zeros(Float64, U.p["N_u"])
    ψ_hor = zeros(Float64, U.p["N_u"])
    ψ_r   = zeros(Float64, U.p["N_u"])
    t_r   = zeros(Float64, U.p["N_u"])
    #--------------------------------------
    #--------------------------------------
    @showprogress for (i,u) in enumerate(U.p["u_grid"])
        U = Update(U)
        ψ_scri[i] = U.ψ[end]
        ind_extract = findall(x -> x > 2.0 * U.p["r_ext"] + u, U.p["v_grid"])[1]
        ψ_r[i] = U.ψ[ind_extract]
        t_r[i] = 0.5 * (u + U.p["v_grid"][ind_extract])
    end
    #--------------------------------------
    #--------------------------------------
    ψ_hor = U.ψ
    #--------------------------------------
    open(filename_H, "w") do f
        write(f, "# v, ψ_hor \n")
        for i = 1:U.p["N_v"]
            write(f, string(U.p["v_grid"][i], "\t", ψ_hor[i], "\n"))
        end
    end
    open(filename_I, "w") do f
        write(f, "# u, ψ_scri \n")
        for i = 1:U.p["N_u"]
            write(f, string(U.p["u_grid"][i], "\t", ψ_scri[i], "\n"))
        end
    end
    open(filename_R, "w") do f
        write(f, "# t, ψ_r \n")
        for i = 1:U.p["N_u"]
            write(f, string(t_r[i], "\t", ψ_r[i], "\n"))
        end
    end
    #--------------------------------------
    #--------------------------------------
    tprint(Panel("Evolution Complete!", box = :ROUNDED, fit = true))
    return nothing
end

end