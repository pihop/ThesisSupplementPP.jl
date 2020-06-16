Directory contains the code for Chapter 6 of the thesis.
* The reaction network model is provided in `rn.jl`. Full stochastic model given 
  by `rn_full`. The dynamics separated into modes is given by `rn_mode_1` and
  `rn_mode_2`.
* Files `hybrid_fluid.jl`, `hybrid_lna.jl` and `hybrid_mc.jl` provide the
  implementations for the approxiamtion methods described in Chapter 6.
* File `opt_goals.jl` implement reward functions on the mean and first moment of
  the trajectories.
* Run scripts `run_script.jl` and `run_script_fixed.jl` calls the
  approximation methods and uses them to calculate the reward functions for
  different model parametrisations.
  `run_script.jl` samples a new set of parameters for each approximation
  methods.
  `run_script_fixed.jl` uses the same set of model parameters for all
  approximations.
  The results are save into `data/`.
* `run_analysis.jl` calculates the accuracy comparison results presented in
  Chapter 6.
* `plot_cost.jl` plots the cost function in the discussion section.
* `run.sh` provides the glue to do everything in a single script.



