This directory provides the code for the running example of Chapter 4 of the thesis.
  * `rn.jl` provides the chemical reaction network models of the running example.
  `rn_full` gives the full stochastic model while `rn_mode_1` and `rn_mode_2`
  give the dynamics separated into two models.
  * `fluid_hybrid_simulation.jl` and `lna_hybrid_simulation.jl` are used to
  construct and simulate the hybrid fluid and linear noise processes.
  * `matlab/run_matlab.sh` clones the github repository of CERENA matlab toolbox
  from `https://github.com/CERENADevelopers/CERENA.git` and uses it to run the 
  model provided in `matlab/running_example/modelDef_running.m`.
  * `run_simulations.sh` runs all simulation. The results are saved in Julia
  langauges `.jld2` format to `data/` directory.
  * `run_plots.sh` runs the plot scripts. The results are saved in `plots/`
  directory.




