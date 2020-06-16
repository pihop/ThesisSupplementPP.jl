__precompile__()
module ThesisSupplementPP

using LinearAlgebra
using SymEngine
using MacroTools: postwalk
using DiffEqBiological
using Combinatorics
#using SymPy
using Logging

include("moment_closure.jl") 
include("normal_closure.jl") 
include("structures.jl") 

export sym_to_expr, get_stoich_prop, gen_diffusion, gen_noise_covar, gen_covar_syms
export make_symbolic_mc

end # module
