#Example: The Quantum Permutation Group on 4 Elements

using Oscar
using QuantumAutomorphismGroups
using QuantumGB


#To get the relations of the quantum permutation group on 4 elements
relat , relat_by_type, u, A = getQuantumPermutationGroup(4);

#You can interreduce the relations
red_relat = interreduce!_noRedTail(relat)

#You could get the obstruction pairs
#red_relat_obstr is made up of tuples of the form (first_poly, second_poly, Spoly_red, Spoly)
red_relat_obstr = get_obstruction_pairs(red_relat);


#This will get us the vector of spolynomials
red_relat_obstr_spoly = [get_obstruction_pairs(red_relat)[i][3] for i in 1:length(red_relat_obstr)];


#This will get us the first spolynomial
red_relat_obstr_spoly[1]
