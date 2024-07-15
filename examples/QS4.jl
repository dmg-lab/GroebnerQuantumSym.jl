#Example: The Quantum Permutation Group on 4 Elements

#= To get this to run

using Pkg
Pkg.add("Oscar#master")
Pkg.add("https://github.com/dmg-lab/QuantumAutomorphismGroups.jl.git")


=#



using Oscar
using QuantumAutomorphismGroups
using QuantumGB


#To get the relations of the quantum permutation group on 4 elements
relat , relat_by_type, u, A = getQuantumPermutationGroup(4);

relat_ = typeof(relat[1])[]
append!(relat_,relat_by_type[:col_sum])
append!(relat_,relat_by_type[:row_sum])


append!(relat_,relat_by_type[:zero_divisor])
ids = relat_by_type[:idempotent]

r, h = normal_form_with_rep(ids[3],relat_)
r

interreduce!_noRedTail(relat_)

println.(relat)

interreduce!_noRedTail(relat)

#=Check Gb
groebner_basis(interreduce!_noRedTail(relat))

#You can interreduce the relations
red_relat = interreduce!_noRedTail(relat)

#You could get the obstruction pairs
#red_relat_obstr is made up of tuples of the form (first_poly, second_poly, Spoly_red, Spoly)
red_relat_obstr = get_obstruction_pairs(red_relat);


#This will get us the vector of spolynomials
red_relat_obstr_spoly = [get_obstruction_pairs(red_relat)[i][3] for i in 1:length(red_relat_obstr)];


#This will get us the first spolynomial
red_relat_obstr_spoly[1]


#Example reps

rscs = typeof(relat[1])[]
append!(rscs,relat_by_type[:col_sum])
append!(rscs,relat_by_type[:row_sum])
append!(rscs,relat_by_type[:zero_divisor])
ids = relat_by_type[:idempotent]
gb1 = groebner_basis(rscs)
id1 = map(x->normal_form(x,gb1),ids)





id1-u[2,3]cs[3]

   

normal_form_with_rep(id1,rcwel)

