#!/bin/bash
#We run the simulations for population level 100 and 5000 simulation trajectories.
julia --project=../../../Project.toml ssa_simulation.jl 100 5000
julia --project=../../../Project.toml fluid_hybrid_simulation.jl 100 5000
julia --project=../../../Project.toml lna_hybrid_simulation.jl 100 5000
cd matlab
sh run_matlab.sh
