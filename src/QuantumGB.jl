module QuantumGB

using Oscar;

greet() = print("Hello World!")

export 
    normal_form_noRedTail,
    normal_form_with_rep,
    interreduce!_noRedTail,
    get_obstruction_pairs,
    getQuantumRelationsByType,
    lm,
    GeneratingSet


include("ReductionAutomata.jl")
include("gbsteps.jl")
include("symmetric_group_gb.jl")
include("interactive_reduce.jl")
#include("pointed_Bases.jl")


end # module QuantumGB
