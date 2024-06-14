module QuantumGB

using Oscar;
using QuantumAutomorphismGroups;

greet() = print("Hello World!")

export 
    normal_form_noRedTail,
    normal_form_with_rep,
    interreduce!_noRedTail,
    get_obstruction_pairs


include("gbsteps.jl")


end # module QuantumGB
