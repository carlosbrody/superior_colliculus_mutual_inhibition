# This script makes the EDF 9b,c,d panels
# This script works in Julia 1.5 (Note the change!)
include("../A_REWRITE/dprime_analysis.jl")

# This makes panels 9b and 9c
scatterplot()

# To make dprime figure, 9d
# This generates example solutions, so the exact figure is variable
# Might require "using JLD2, FileIO"
save_solutions()
plot_dprime()

