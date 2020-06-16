This is a supplement to my PhD thesis titled "Formal modelling and
approximation-based analysis for mode-switching population dynamics".

It contains the simulation examples used in the thesis.

# Requirements
Julia 1.4 is needed to run the code. Available from https://julialang.org/downloads/.

# Instructions
Easiest way to run the code is to use Julia's built-in package manager to clone
the repository.
For that, in Julia's REPL do
```
pkg> develop https://github.com/pihop/ThesisSupplement.jl
```
This clones this repository to the `dev` directory of the local Julia
installation.

For running example in Chapter 4 see:
```
examples/thesis_ch_4/running_ex/
```
The script `run_simulation.sh` generates the simulation data and `run_plots.sh`
plots the results.

For running example in Chapter 5 see:
```
examples/thesis_ch_5/running_ex/
```
The script `run.sh` generates the simulation data and plots the results.
The results in the results of Chapter 5 are based on examples in
```
examples/thesis_ch_5/maze_ex/
examples/thesis_ch_5/three_mode_ex/
```
These again come with their own run scrips to generate the data.
```
examples/thesis_ch_5/run_long.sh
```
calls both of these.
Finally, 
```
examples/thesis_ch_5/analyse.sh
```
produces the plots appearing in the thesis and the data for the tables (written
to `examples/thesis_ch_5/output.txt`.

For running example in Chapter 6 see:
```
examples/thesis_ch_6/running_ex/
```
The script `run.sh` generates the simulation data and generates the data for the
table.

