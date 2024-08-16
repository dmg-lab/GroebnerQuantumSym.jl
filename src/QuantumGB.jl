module QuantumGB

using Oscar;
using QuantumAutomorphismGroups;

greet() = print("Hello World!")

export 
    normal_form_noRedTail,
    normal_form_with_rep,
    interreduce!_noRedTail,
    get_obstruction_pairs,
    getQuantumRelationsByType,
    lm,
    GeneratingSet


include("gbsteps.jl")
include("GeneratingSet.jl")
#include("pointed_Bases.jl")


end # module QuantumGB
