# Julia script for aggregating the analyis of the results.
# CONSTANTS
# number of parametrisations the simulations were run for.
iters = 50 


include("results.jl")
include("plots.jl")

# Table 5.2 Timings for the maze example
println("Timings for the maze example")
all_timings("maze_ex", iters)

# Table 5.3
println("Mean maximum error analysis (for mean statistics)")
println("MAZE EXAMPLE")
# idx 13 corresponds to the location (3,3) in the model.
all_errors_means("maze_ex", 13, iters)
println()
println("THREE MODE EXAMPLE")
# idx 4 corresponds to the location (1,1) in the model.
all_errors_means("three_mode_ex", 3, iters)
println()

# Table 5.3
println("Mean maximum error analysis (for variance statistics)")
println("MAZE EXAMPLE")
# idx 13 corresponds to the location (3,3) in the model.
all_errors_vars("maze_ex", 13, "lna", iters)
println()
println("THREE MODE EXAMPLE")
# idx 3 corresponds to the location (1,1) in the model.
all_errors_vars("three_mode_ex", 3, "lna", iters)
println()

# Generate all of the plots from Chapter 5.
# Plots are written to (three_mode_ex or maze_ex)/plots/
# Figures 5.10
make_heat_plot("three_mode_ex", 100, 3, iters)
make_heat_plot("three_mode_ex", 200, 3, iters)

# Figures 5.11
make_heat_plot("maze_ex", 100, 13, iters)
make_heat_plot("maze_ex", 200, 13, iters)

# Figures 5.12 and 5.14. The plots are writtent to (three_mode_ex or
# maze_ex)/plots/var_100_i.pdf
# The corresponding parameter values for the model are in (three_mode_ex or
# maze_ex)/plots/var_100_i.txt

for i in 1:50
    make_var_plot("three_mode_ex", 100, i, [(1,2)], "lna", iters, ylim=(0.0, 0.2))
    make_var_plot("maze_ex", 100, i, [(1,4)], "lna", iters, ylim=(0.0, 0.3))
end

# Figures 5.13
make_heat_plot_var("three_mode_ex", 100, 3, "lna", iters)
make_heat_plot_var("three_mode_ex", 200, 3, "lna", iters)

# Figures 5.15
make_heat_plot_var("maze_ex", 100, 13, "lna", iters)
make_heat_plot_var("maze_ex", 200, 13, "lna", iters)


