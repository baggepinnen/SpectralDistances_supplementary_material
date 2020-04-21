# Supplamentary material to SpectralDistances.jl
This repository contains the code to generate the figures in the paper

["New Metrics Between Rational Spectra and their Connection to Optimal Transport", Bagge Carlson and Chitre 2020](http://arxiv.org/abs/2004.09152)

with the accompanying software package [SpectralDistances.jl](https://github.com/baggepinnen/SpectralDistances.jl).

The code is written in Julia, https://julialang.org/

## Installation instructions
1. Install julia v1.4 or later ([instructions](https://julialang.org/downloads/))
2. To ensure that the same environment is instantiated as was used to create the figures in the paper, use the Manifest.toml file checked into this repository. Natigate to the folder of this repo, start julia, then run
```julia
julia> using Pkg
julia> pkg"activate ."
julia> pkg"instantiate -m"
```
3. Examples can now be run by including the corresponding script. Example:
```julia
julia> include("kbarycenters_birds.jl")
```




## Datasets
The code assumes that the datasets provided under the links below exist in subfolders to this folder, named
- `./birds`
- `./ships`

### Links
- [Birds data: 2.5GB train, 356 MB test. (google drive)](https://drive.google.com/open?id=1jMGMjj-KPJ8b4qoFo5WZqrOVr3sHzwDf)

## Scripts
- `cumulative_spectra.jl` reproduces figure 2.
- `varying_parameters.jl` reproduces figures 4-5.
- `interpolation.jl` reproduces figures 8-9.
- `ship_detection.jl` reproduces figures 10-11.
- `kbarycenters_birds.jl` reproduces figure 12.
