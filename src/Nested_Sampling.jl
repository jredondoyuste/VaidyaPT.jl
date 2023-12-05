module Nested_Sampling

using NestedSamplers
using Distributions
using LinearAlgebra
using StatsBase
using DelimitedFiles
using Parameters
using Interpolations
using Noise

export SimData, Inference_Config, inference, read_data, shift, cut, shift_cut, reconstruct, get_vshift, add_noise

"""
    SimData(v, ψ): Simulated data structure.
"""
struct SimData
    v::Vector{Float64}
    ψ::Vector{Float64}
end

function SimData(parent::String, m2::Float64, τ::Float64)
    read_data(parent * "data/m2_$(m2)_tau_$(τ)")
end

"""
    read_data(dir): Read the data from a file and return a SimData struct.
"""
function read_data(dir)
    fn = joinpath(dir, "Horizon_l_2_s_2.dat")
    v = readdlm(fn, '\t', Float64, comments = true, comment_char = '#')[:, 1]
    ψ = readdlm(fn, '\t', Float64, comments = true, comment_char = '#')[:, 2]
    SimData(v, ψ)
end

"""
    add_noise(s, σ): Add Gaussian noise to a waveform.
"""
function add_noise(s::SimData, σ::Float64)
    ψ = add_gauss(s.ψ, σ)
    SimData(s.v, ψ)
end

"""
    get_vshift(s): Get the time shift of the waveform so that v = 0 corresponds to the peak.
"""
function get_vshift(s::SimData)
    s.v[argmax(abs.(s.ψ))]
end

"""
    shift(s): Shift the waveform so that v = 0 corresponds to the peak.
"""
function shift(s::SimData)
    vpeak = s.v[argmax(abs.(s.ψ))]
    v = s.v .- vpeak
    SimData(v, s.ψ)
end

"""
    cut(s, v0, vend): Cut the waveform between v0 and vend.
"""
function cut(s::SimData, v0::Float64, vend::Float64)
    i0 = findfirst(s.v .> v0)
    iend = findfirst(s.v .> vend)
    SimData(s.v[i0:iend], s.ψ[i0:iend])
end

"""
    shift_cut(s, v0, vend): Shift and cut the waveform between v0 and vend.
"""
function shift_cut(s::SimData, v0::Float64, vend::Float64)
    cut(shift(s), v0, vend)
end

"""
    get_error(s_high, s_low): Get the error between two waveforms evaluated at different resolutions.
"""
function get_error(s_high::SimData, s_low::SimData)
    ψ2 = linear_interpolation(s_low.v, s_low.ψ, extrapolation_bc=Line())(s_high.v)
    return abs.(s_high.ψ - ψ2)
end

"""
    Inference_Config(data_dir, v0, vend, err_type, err_constant, lower_res, nlive, model, m, vchange, α, β, tau, ln_A_min, ln_A_max, ω_min, ω_max, φ_min, φ_max, τ_min, τ_max): Configuration for the nested sampling inference.
"""
@with_kw struct Inference_Config
    data_dir::String
    v0::Float64 = 10.0
    vend::Float64 = 250.0
    err_type::String = "constant"
    err_constant::Float64 = 1e-3
    lower_res::String = data_dir * "_low"
    nlive::Int64 = 16
    model::String = "ds"
    m::Function = x -> 1.0
    vchange::Float64 = 100.0
    α::Float64 = 4.19
    β::Float64 = 5.03
    tau::Float64 = 10.0
    ln_A_min::Float64 = -10.0
    ln_A_max::Float64 = 2.0
    ω_min::Float64 = 0.0
    ω_max::Float64 = 1.0
    φ_min::Float64 = 0.0
    φ_max::Float64 = 2π
    τ_min::Float64 = 1.0
    τ_max::Float64 = 25.0
end

"""
    oneds(t, pars): One damped sinusoid waveform
"""
function oneds(t, pars)
    lnA = pars[1]
    ω = pars[2]
    φ = pars[3]
    τ = pars[4]
    return @. exp(lnA) * exp( - t / τ) * cos(ω * t + φ)
end

"""
    dyn_RD(t, pars, m, v1, α, β): Dynamical ringdown waveform
"""
function dyn_RD(t, pars, m, v1, α, β)
    lnA = pars[1]
    φ = pars[2]
    f_qnm = 0.373672
    tau_qnm = 11.2407
    m_end = m(1000)
    m_0 = m(0)
    ω0 = f_qnm / m_0 - im / (tau_qnm * m_0)
    ωend = f_qnm / m_end - im / (tau_qnm * m_end)
    dm = m_end - m_0
    q_factor = (1 + α * dm * exp(  im * β * dm)) * exp(im * ωend * v1) * exp(-im * ω0 * v1) 
    ratio_amp = abs(q_factor)
    delta_phase = angle(q_factor)
    return @. exp(lnA) * ( 1 + (ratio_amp -1) * (m(t) - m_0) / dm) * exp( - t / (tau_qnm * m(t))) * cos(f_qnm * t / m(t) + φ - delta_phase * (m(t) - m_0) / dm)
end

"""
    likelihood_ds(t, ψ, err_type, err_level, err_ψ): Likelihood function for the damped sinusoid model.
"""
function likelihood_ds(t, ψ, err_type, err_level, err_ψ)
    lh_c(X) = - 0.5 * sum( ((ψ - oneds(t, X)) ./ (1e-16 + err_level)).^2)
    lh_r(X) = - 0.5 * sum( ((ψ - oneds(t, X)) ./ (1e-16 .+ err_ψ)).^2)
    if err_type == "constant"
        return lh_c
    elseif err_type == "resolution"
        return lh_r
    elseif err_type == "noise"
        return lh_c
    end
end

"""
    likelihood_dr(t, ψ, m, v1, α, β, err_type, err_level, err_ψ): Likelihood function for the dynamical ringdown model.
"""
function likelihood_dr(t, ψ, m, v1, α, β, err_type, err_level, err_ψ)
    lh_c(X) = - 0.5 * sum( ((ψ - dyn_RD(t, X, m, v1, α, β)) ./ (1e-16 + err_level)).^2)
    lh_r(X) = - 0.5 * sum( ((ψ - dyn_RD(t, X, m, v1, α, β)) ./ (1e-16 .+ err_ψ)).^2)
    if err_type == "constant"
        return lh_c
    elseif err_type == "resolution"
        return lh_r
    elseif err_type == "noise"
        return lh_c
    end
end

"""
    inference(i): Perform the nested sampling inference.
"""
function inference(i::Inference_Config)
    sd = read_data(i.data_dir)
    sd_low = SimData([0.0], [0.0])
    sd = shift_cut(sd, i.v0, i.vend) 
    err = zeros(length(sd.ψ))
    if i.err_type == "resolution"
        sd_low = read_data(i.lower_res)
        sd_low = shift_cut(sd_low, i.v0, i.vend)
        err = get_error(sd, sd_low)
    elseif i.err_type == "noise"
        sd = add_noise(sd, i.err_constant)
    end
    if i.model == "ds"
        priors = [
            Uniform(i.ln_A_min, i.ln_A_max),
            Uniform(i.ω_min, i.ω_max),
            Uniform(i.φ_min, i.φ_max),
            Uniform(i.τ_min, i.τ_max)
        ]
        lh = likelihood_ds(sd.v, sd.ψ, i.err_type, i.err_constant, err)
        model = NestedModel(lh, priors)
        sampler = Nested(4, i.nlive)
        chain, state = sample(model, sampler; dlogz=0.5, param_names = ["lnA", "ω", "φ", "τ"])
        chain = sample(chain, Weights(vec(chain["weights"])), length(chain))
        return chain, state
    elseif i.model == "dr"
        priors = [
            Uniform(i.ln_A_min, i.ln_A_max),
            Uniform(i.φ_min, i.φ_max)
        ]
        lh = likelihood_dr(sd.v, sd.ψ, i.m, i.vchange, i.α, i.β, i.err_type, i.err_constant, err)
        model = NestedModel(lh, priors)
        sampler = Nested(2, i.nlive)
        chain, state = sample(model, sampler; dlogz=0.5, param_names = ["lnA", "φ"])
        println("\n")
        chain = sample(chain, Weights(vec(chain["weights"])), length(chain))
        println("\n")
        return chain, state
    else
        println("Model not supported")
    end
end

"""
    reconstruct(v, chain, model; m, v1, α, β): Reconstruct the waveform from the chain.
"""
function reconstruct(v, chain, model; m = x -> 0.0, v1 = 100.0, α = 4.19, β = 5.03)
    if model == "ds"
        ψ_reconstructed =  oneds(v, [mean(chain[:lnA]), mean(chain[:ω]), mean(chain[:φ]), mean(chain[:τ])])
        quantiles = [quantile(collect(chain[:lnA][:,1]), [0.95, 0.05]), quantile(collect(chain[:ω][:,1]), [0.95, 0.05]), quantile(collect(chain[:φ][:,1]), [0.95, 0.05]), quantile(collect(chain[:τ][:,1]), [0.95, 0.05])]
        ψ_5 = oneds(v, [quantiles[1][1], quantiles[2][1], quantiles[3][1], quantiles[4][1]])
        ψ_95 = oneds(v, [quantiles[1][2], quantiles[2][2], quantiles[3][2], quantiles[4][2]])
        return ψ_reconstructed, ψ_5, ψ_95
    elseif model == "dr"
        ψ_reconstructed =  dyn_RD(v, [mean(chain[:lnA]), mean(chain[:φ])], m, v1, α, β)
        quantiles = [quantile(collect(chain[:lnA][:,1]), [0.95, 0.05]), quantile(collect(chain[:φ][:,1]), [0.95, 0.05])]
        ψ_5 = dyn_RD(v, [quantiles[1][1], quantiles[2][1]], m, v1, α, β)
        ψ_95 = dyn_RD(v, [quantiles[1][2], quantiles[2][2]], m, v1, α, β)
        return ψ_reconstructed, ψ_5, ψ_95
    end
end

end