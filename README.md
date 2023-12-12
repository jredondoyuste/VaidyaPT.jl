# VaidyaPT

[![Build Status](https://github.com/jredondoyuste/VaidyaPT.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jredondoyuste/VaidyaPT.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/jredondoyuste/VaidyaPT.jl.svg?branch=main)](https://travis-ci.com/jredondoyuste/VaidyaPT.jl)
[![Coverage](https://codecov.io/gh/jredondoyuste/VaidyaPT.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jredondoyuste/VaidyaPT.jl)

This is an auxiliary repository containing the source code used in [arXiv:2312.XXXXX](URL). 

The code includes scripts to evolve axial gravitational perturbations in the Vaidya spacetime, using a characteristic algorithm.
Moreover, utilities to analyze the data are included: both an algorithm to minimize the mismatch between the waveform an a model consisting on superposed damped sinusoids with constant amplitudes and frequencies, as well as a Bayesian-based inference routine, which can be used to compare a model based on damped sinusoids, with a _dynamical ringdown_ model.

Most of the functionalities of the code are shown in the [`example.jl`](example/example.jl) script. 

## Installation and usage

1.-Install [Julia](https://julialang.org/downloads/) 
2.-Clone this repository:
```
git clone https://github.com/jredondoyuste/VaidyaPT.jl.git
```
3.-Launch julia and install the code. You will also need the `Plots.jl` package to run the example script (but this can be easily replaced by your favourite plotting package).
```
julia --project=VaidyaPT.jl
using Pkg;
Pkg.add("Plots");
Pkg.instantiate();
Pkg.precompile();
```

Then, you should be able to run the code as you wish. To get started, check that the [`example.jl`](example/example.jl) script runs:

```
julia --project=VaidyaPT.jl VaidyaPT/example/example.jl
```

This should take a few minutes to run. After it finishes, you can see the results in the `example/figures/` directory.

## Questions and Contributions

The code can be seeen as a scheme where the master equations in more complicated spacetimes can be added with some effort. If you are interested in extending the code, or find a bug, please submit a PR or get in touch at jaime.redondo.yuste[at]nbi[dot]ku[dot]dk.

## Citing

If you make use of this code, please consider citing the companion paper using the [`CITATION.bib`](CITATION.bib) template:

```
@article{Redondo-Yuste:2023ipg,
    author = "Redondo-Yuste, Jaime and Pere\~niguez, David and Cardoso, Vitor",
    title = "{Ringdown of a dynamical spacetime}",
    eprint = "2312.04633",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "12",
    year = "2023"
}
```
