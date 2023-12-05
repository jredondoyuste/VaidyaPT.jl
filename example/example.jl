using VaidyaPT              # import the package 
using Plots                 # import the plotting backend
using LaTeXStrings          # import LaTeXStrings for LaTeX labels
using StatsBase             # import StatsBase for the mean function

m2_val = 1.05
if !isdir("example/data/m2_1.05")
    mkdir("example/data/m2_1.05")
end

Ψ = Field(parfile = "example/Schwarzschild_pars.txt")     # create a field object
Initialize(Ψ)                                             # initialize the field 
Evolve(Ψ, "example/data/Schwarzschild/")                  # evolve the field and save the data in the folder "examples/data/"
p = read_pars("example/Schwarzschild_pars.txt")             # read the parameters
p = change_par(p, "m_2", 1.05)                              # change the final mass of the black hole
Ψ = Field(p = p)                                            # create a new field object with the new parameters
Initialize(Ψ)                                               # initialize the field
Evolve(Ψ, "example/data/m2_1.05/")                          # evolve the field and save the data in the folder "examples/data/"

# plot the data

data_S = read_data("example/data/Schwarzschild");
data_V = read_data("example/data/m2_1.05");

plot(data_S.v, abs.(data_S.ψ) .+ 1e-30, label = L"m_2 = m_1", xlabel = L"v", ylabel = L"\Psi(ℋ)", legend = :topleft)
plot!(data_V.v, abs.(data_V.ψ) .+ 1e-30, label = L"m_2 = 1.05 m_1", xlabel = L"v", ylabel = L"\Psi(ℋ)", legend = :topleft)
plot!(yscale = :log10)
ylims!(1e-12,1)
xlims!(0.0, 150.0)
savefig("example/Evolution.pdf")

# fit the data with damped sinusoids 

ic_1 = Inference_Config(
	data_dir = "example/data/m2_1.05",
	v0 = 10.0,
	vend = 150.0,
	err_type = "noise",
	err_constant = 1e-4,
	nlive = 16,
)

chain1, state1 = inference(ic_1)

# fit the data with the dynamical ringdown model

m2 = 1.05
tau = 10.0
v1 = 100.0 - get_vshift(data_V)
α_factor = 4.19
β_factor = 5.03

function mass_fun(x)
    return 1.0 + 0.5 * (m2 - 1.0) * (1 + tanh((x - v1) / tau))
end

ic_2 = Inference_Config(
	data_dir = "example/data/m2_1.05",
	v0 = 10.0,
	vend = 150.0,
	err_type = "noise",
	err_constant = 1e-4,
	model = "dr",
	m = mass_fun,
	vchange = v1,
	α = α_factor,
	β = β_factor,
	tau = 10.0,
	nlive = 16,
)

chain2, state2 = inference(ic_2)

# plot the reconstructed waveforms

data_cut = shift_cut(data_V, 10.0, 150.0)
data_noisy = add_noise(data_cut, 1e-4)
ψ1, _, _ = reconstruct(data_noisy.v, chain1, "ds")
ψ2, _, _ = reconstruct(data_noisy.v, chain2, "dr"; m = mass_fun, v1 = v1, α = α_factor, β = β_factor)

plot(data_noisy.v, abs.(data_noisy.ψ) .+ 1e-30, label = L"Injection", xlabel = L"v", ylabel = L"\Psi(ℋ)", legend = :topright, color = :gray)
plot!(data_cut.v, abs.(data_cut.ψ) .+ 1e-30, label = L"Data", color= :black)
plot!(data_noisy.v, abs.(ψ1) .+ 1e-30, label = L"DS", color = :blue)
plot!(data_noisy.v, abs.(ψ2) .+ 1e-30, label = L"DR", color = :red)

plot!(yscale = :log10)
ylims!(1e-7,0.2)
xlims!(0.0, 150.0)

savefig("example/Reconstruction.pdf")