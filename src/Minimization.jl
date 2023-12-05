module Minimization

export mismatch, minimize, find_N_cycles

using DelimitedFiles
using LineSearches
using Optim

"""
    wf_ds(t, pars): Damped sinusoid waveform.
"""
function wf_ds(t, pars)
    Amp = pars[1]
    freq = pars[2]
    phase = pars[3]
    damping = pars[4]
    return @. Amp * exp(-t / damping) * cos(freq * t + phase)
end

"""
    wf_2ds(t, pars): Sum of two damped sinusoids.
"""
function wf_2ds(t, pars)
    Amp = pars[1]
    freq = pars[2]
    phase = pars[3]
    damping = pars[4]
    Amp2 = pars[5]
    freq2 = pars[6]
    phase2 = pars[7]
    damping2 = pars[8]
    return @. Amp * exp(-t / damping) * cos(freq * t + phase) +
              Amp2 * exp(-t / damping2) * cos(freq2 * t + phase2)
end

"""
    wf_qnm_ds(t, pars, f0, tau0): Sum of two damped sinusoids, with the frequencies of the first fixed.
"""
function wf_qnm_ds(t, pars, f0, tau0)
    Amp0 = pars[1]
    phase0 = pars[2]
    Amp1 = pars[3]
    phase1 = pars[4]
    f1 = pars[5]
    tau1 = pars[6]
    return @. Amp0 * exp(-t / tau0) * cos(f0 * t + phase0) +
              Amp1 * exp(-t / tau1) * cos(f1 * t + phase1)
end

"""
    inner_product(t, psi1, psi2): Inner product of two waveforms.
"""
function inner_product(
    t::Vector{Float64},
    psi1::Vector{Float64},
    psi2::Vector{Float64},
)::Float64
    return sum(psi1 .* psi2) * (t[2] - t[1])
end

"""
    mismatch(t, psi1, psi2): Mismatch between two waveforms.
"""
function mismatch(
    t::Vector{Float64},
    psi_1::Vector{Float64},
    psi_2::Vector{Float64},
)::Float64
    return 1 -
           inner_product(t, psi_1, psi_2) /
           sqrt(inner_product(t, psi_1, psi_1) * inner_product(t, psi_2, psi_2))
end

wf_models = Dict(
    "ds" => wf_ds,
    "2ds" => wf_2ds,
    "qnm_ds" => wf_qnm_ds,
)

pars_models = Dict(
    "ds" => ["A", "f", "phi", "τ"],
    "2ds" => ["A0", "f0", "phi0", "τ0", "A1", "f1", "phi1", "τ1"],
    "qnm_ds" => ["A0", "phi0", "A1", "phi1", "f1", "τ1"],
)

"""
    minimize(wf_model, data_file, t0, tend; f0, tau0, tlim, tol, max_iters, Amp_guess): Minimize the mismatch between a waveform model and a data file.
"""
function minimize(
    wf_model::String,
    data_file::String,
    t0::Float64,
    tend::Float64;
    f0::Float64 = 0.373672,
    tau0::Float64 = 11.2407,
    tlim::Float64 = 10.0,
    tol::Float64 = 1e-3,
    max_iters::Int = 10,
    Amp_guess::Float64 = 0.1,
)
    # --------------------------------------
    mism_val = 1.0
    iters = 0
    # --------------------------------------
    data = readdlm(data_file, '\t', Float64, comments = true, comment_char = '#')
    t_full = data[:, 1]
    psi_full = data[:, 2]
    tshift = t_full[argmax(abs.(psi_full))]
    t_shifted = t_full .- tshift
    t = t_shifted[(t_shifted.>=t0).&(t_shifted.<=tend)]
    psi = psi_full[(t_shifted.>=t0).&(t_shifted.<=tend)]
    # --------------------------------------
    wf_fun = wf_models[wf_model]
    wf_ds_fun(x, betas) = wf_fun(x, betas)
    wf_2ds_fun(x, betas) = wf_fun(x, betas)
    wf_qnm_ds_fun(x, betas) = wf_fun(x, betas, f0, tau0)
    if wf_model == "ds"
        wf = wf_ds_fun
    elseif wf_model == "2ds"
        wf = wf_2ds_fun
    elseif wf_model == "qnm_ds"
        wf = wf_qnm_ds_fun
    else
        error("Waveform model $(wf_model) not implemented")
    end
    model_pars = pars_models[wf_model]
    # --------------------------------------
    function sqerror(betas, X, Y)
        err = 0.0
        for (i, x) in enumerate(X)
            pred_i = wf(x, betas)
            err += (Y[i] - pred_i)^2
        end
        return err
    end
    # --------------------------------------
    while (mism_val > tol) 
        if iters > max_iters
            break
        end
        iters += 1
        # --------------------------------------
        ini_guess = zeros(Float64, length(model_pars))
        for (i, pp) in enumerate(model_pars)
            if contains(pp, "A")
                ini_guess[i] = Amp_guess + randn() * 0.1 * Amp_guess
            elseif contains(pp, "phi")
                ini_guess[i] = randn() * 0.1
            elseif contains(pp, "f")
                ini_guess[i] = 0.5 + randn() * 0.1 * f0
            elseif contains(pp, "τ")
                ini_guess[i] = 10.0 + randn() * 0.1 * tau0
            end
        end
        # --------------------------------------
        fn = optimize(
            b -> sqerror(b, t, psi),
            ini_guess,
            LBFGS(linesearch = LineSearches.BackTracking()),
            Optim.Options(
                x_tol = 1e-30,
                f_tol = 1e-30,
                g_tol = 1e-30,
                iterations = 100000,
                time_limit = tlim,
            ),
            autodiff = :forward,
        )
        global best_fit = Optim.minimizer(fn)
        global reconstructed = wf(t, best_fit)
        mism_val = mismatch(t, psi, reconstructed)
    end
    return best_fit, reconstructed, mism_val
end

"""
    find_N_cycles(data_file, t0, Ncycles): Find the time at which the Ncycles-th cycle of the waveform starts.
"""
function find_N_cycles(data_file::String, t0::Float64, Ncycles::Int64)::Float64
    data = readdlm(data_file, '\t', Float64, comments = true, comment_char = '#')
    t_full = data[:, 1]
    psi_full = data[:, 2]
    tshift = t_full[argmax(abs.(psi_full))]
    t_shifted = t_full .- tshift
    t = t_shifted[(t_shifted.>=t0)]
    psi = psi_full[(t_shifted.>=t0)]
    tmax = t[1]
    crossings = 0
    for (i, t_i) in enumerate(t)
        if psi[i] * psi[i+1] < 0
            crossings += 1
            if crossings == Ncycles
                tmax = t_i
                break
            end
        end
    end
    return tmax
end

end