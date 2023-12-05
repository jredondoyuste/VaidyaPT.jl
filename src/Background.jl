module Background

export mass!, r_sol!, f_sol!

using Roots
using DifferentialEquations

const half = 1.0 / 2.0
const third = 1.0 / 3.0
const twelfth = 1.0 / 12.0
const r_switch_method = 5.0 / 2.0

"""
    mass!(v, m1, m2, v1, τ): Mass profile of the accreted matter. The accretion interpolates between masses m1 and m2, occurs centered at a time v1, for a timescale τ.
"""
function mass!(
    v::Union{Float64,Vector{Float64}},
    m1::Float64,
    m2::Float64,
    v1::Float64,
    τ::Float64,
    )::Union{Float64,Array{Float64}}
    return @. m1 + half * (m2 - m1) * (1 + tanh((v - v1) / τ))
end

"""
    radial_eq!(r, m): Right hand side of the radial transport equation.
"""
function radial_eq!(r, m)
    return half - m / r
end

"""
    r_ic!(u0, v0, m, tol): Implementation of the initial condition for the radial transport equation (inverting the equation for the tortoise coordinate)
"""
function r_ic!(
    u0::Float64, 
    v0::Float64, 
    m::Float64, 
    tol::Float64 = 1e-13,
    )::Float64

    f(x) = x + 2 * m * log(abs(x / (2 * m) - 1)) - half * (v0 - u0) 
    r_star = half * (v0 - u0)
    if r_star > r_switch_method * m
        guess = r_star
        zero = find_zero(f, guess)
    else
        guess = r_star
        delta = 1
        while delta > tol
            r_new = 2 * m * (1 + exp((r_star - guess) / (2 * m)))
            delta = abs(r_new - guess)
            guess = r_new
        end
        if guess == 2.0 * m
            zero = 2.0 * m + 1e-2 * tol
        else
            zero = guess
        end
    end
    return zero
end

"""
    r_sol!(u0, m1, m2, v0, v1, vend, τ, Δv): Solution to the radial transport equation, i.e., profile of r(u=u0,v), for the mass profile fixed by the input parameters, where the initial conditions are chosen teleologically at vend.  
"""
function r_sol!(
    u0::Float64,
    m1::Float64,
    m2::Float64,
    v0::Float64,
    v1::Float64,
    vend::Float64,
    τ::Float64,
    Δv::Float64,
    )::Vector{Float64}
    ode_problem(u, p, t) = radial_eq!(u, mass!(t, m1, m2, v1, τ))
    ini = r_ic!(u0, vend, m2)
    prob = ODEProblem(ode_problem, ini, (vend, v0), saveat = Δv)
    sol = solve(prob, RadauIIA5(), abstol = 1e-20 * Δv, reltol = 1e-20 * Δv)
    return reverse(sol.u)
end

"""
    f_sol!(u0, m1, m2, v0, v1, vend, τ, Δv, Δu): Obtain the background function f(u=u0, v) by solving the background Einstein equation.
"""
function f_sol!(
    u0::Float64,
    m1::Float64,
    m2::Float64,
    v0::Float64,
    v1::Float64,
    vend::Float64,
    τ::Float64,
    Δv::Float64,
    Δu::Float64,
    )::Vector{Float64}
    r_2 = r_sol!(u0 + 2.0 * Δu, m1, m2, v0, v1, vend, τ, Δv)
    r_m2 = r_sol!(u0 - 2.0 * Δu, m1, m2, v0, v1, vend, τ, Δv)
    r_1 = r_sol!(u0 + Δu, m1, m2, v0, v1, vend, τ, Δv)
    r_m1 = r_sol!(u0 - Δu, m1, m2, v0, v1, vend, τ, Δv)
    return @. -(twelfth * r_m2 - 2.0 * third * r_m1 + 2.0 * third * r_1 - twelfth * r_2) /
              Δu
end

end