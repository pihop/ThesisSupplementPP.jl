#!/bin/bash
#julia --project=../../Project.toml run_script.jl 100
julia --project=../../Project.toml run_script_fixed.jl 100
julia --project=../../Project.toml run_script_fixed.jl 200
#julia --project=../../Project.toml run_analysis.jl "data/data_100.jld2" > output.txt
julia --project=../../Project.toml run_analysis.jl "data/data_f100.jld2" > output_f100.txt
#julia --project=../../Project.toml run_analysis.jl "data/data_200.jld2" > output.txt
julia --project=../../Project.toml run_analysis.jl "data/data_f200.jld2" > output_f200.txt
#julia --project=../../Project.toml plot_cost.jl 
