module PerturbationEquations

export Potential!, Source!

"""
    Source!(r, f, ϕ): Compute the source term of the perturbation equation.
"""
function Source!(
    r::Union{Float64,Array{Float64}},
    f::Union{Float64,Array{Float64}},
    ϕ::Union{Float64,Array{Float64}},
)
    return @. - f * ϕ / (r^2)
end

"""
    Potential!(s, l, m, r): Compute the potential term of the perturbation equation.
"""
function Potential!(
    s::Int64,
    l::Int64,
    m::Union{Float64,Array{Float64}},
    r::Union{Float64,Array{Float64}},
)
    return @. l * (l + 1) / (2 * r^2) + (1 - s^2) * m / r^3
end

end