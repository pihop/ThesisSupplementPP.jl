#!/bin/bash
# Run script gathering the simulation and plotting.
julia --project=../../../Project.toml hybrid_fluid.jl
julia --project=../../../Project.toml hybrid_lna.jl
julia --project=../../../Project.toml ssa_inhomogeneous.jl
julia --project=../../../Project.toml ssa.jl 5000
julia --project=../../../Project.toml plots.jl
