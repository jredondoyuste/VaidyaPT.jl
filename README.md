# VaidyaPT

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jredondoyuste.github.io/VaidyaPT.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jredondoyuste.github.io/VaidyaPT.jl/dev/)
[![Build Status](https://github.com/jredondoyuste/VaidyaPT.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jredondoyuste/VaidyaPT.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/jredondoyuste/VaidyaPT.jl.svg?branch=main)](https://travis-ci.com/jredondoyuste/VaidyaPT.jl)
[![Coverage](https://codecov.io/gh/jredondoyuste/VaidyaPT.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jredondoyuste/VaidyaPT.jl)

This is an auxiliary repository containing the source code used in [arXiv:2312.XXXXX](URL). 

The code includes scripts to evolve axial gravitational perturbations in the Vaidya spacetime, using a characteristic algorithm.
Moreover, utilities to analyze the data are included: both an algorithm to minimize the mismatch between the waveform an a model consisting on superposed damped sinusoids with constant amplitudes and frequencies, as well as a Bayesian-based inference routine, which can be used to compare a model based on damped sinusoids, with a _dynamical ringdown_ model.

Most of the functionalities of the code are shown in the [`example.jl`](example/example.jl) script. 

## Installation and usage

Install Julia (the code has only been tested in Julia 1.8.5) and run

```
]add https://github.com/jredondoyuste/VaidyaPT.jl.git
```

Now, check that you can run the [`example.jl`](example/example.jl) script. The first time (due to compilation time) it might be somewhat slower.

## Questions and Contributions

The code can be seeen as a scheme where the master equations in more complicated spacetimes can be added with some effort. If you are interested in extending the code, or find a bug, please submit a PR or get in touch at jaime.redondo.yuste[at]nbi[dot]ku[dot]dk.

## Citing

If you make use of this code, please consider citing the companion paper using the [`CITATION.bib`](CITATION.bib) template:

```
@article{Redondo-Yuste:2023abc,
    author = "Redondo-Yuste, Jaime and Pereniguez, David and Cardoso, Vitor",
    title = "{Ringdown on a dynamical spacetime}",
    eprint = "2312.XXXXX",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "",
    year = "2023"
}
```
