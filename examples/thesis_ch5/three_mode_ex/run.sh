#!/bin/bash
#julia --project=../../../Project.toml mc_eqns.jl 

# parameter: pop size, simulation end time, number of ssa samples, number of random parameters.
julia --project=../../../Project.toml run_script.jl 100 100.0 5000 50
julia --project=../../../Project.toml run_script.jl 200 100.0 5000 50
julia --project=../../../Project.toml run_script.jl 300 100.0 5000 50 
julia --project=../../../Project.toml run_script.jl 500 10.0 5000 50 
julia --project=../../../Project.toml run_script.jl 1000 10.0 5000 50
