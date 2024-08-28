using Base: GenericIOBuffer
using Oscar: leading_monomial
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




n = 4

rel, relat_by_type, u, A = getQuantumPermutationGroup(n)
cs = relat_by_type[:col_sum]
rs = relat_by_type[:row_sum][2:end]

ip = relat_by_type[:idempotent]
nonip = filter(x->leading_monomial(x) in vcat([u[1,j]^2 for j in 1:n],[u[i,1]^2 for i in 1:n]),ip)
ip_f = filter(x->!(x in nonip),ip)

wels = [u[i,j]*u[i,k] for i in 2:n, j in 2:n, k in 2:n if j != k]
inj = [u[i,j]*u[k,j] for i in 2:n, j in 2:n, k in 2:n if i != k]

rwels = reshape([rwel(k,j,n) for j in 2:n, k in 2:n if j!=k],(n-1)*(n-2))
rinjs = reshape([rinj(k,j,n) for j in 2:n, k in 2:n if j!=k],(n-1)*(n-2))

#filter out shit
rwels_min = filter(x->lm(x) != u[2,2]*u[3,3],rwels)

bg1 = filter(x-> lm(x) ==u[2,2]*u[3,3],rwels)[1]

gens = vcat(ip_f,rs,cs,wels,inj,rwels_min,rinjs)

normal_form_noRedTail(bg1,gens)



gb1 = groebner_basis(gens)

Oscar.AbstractAlgebra.interreduce!(gb1)

gb2 = groebner_basis(rel)

Oscar.AbstractAlgebra.interreduce!(gb2)

gb2[end] in gb1

map(x -> x in gb2,gb1)
map(y -> y in gb1, gb2)

sum(normal_form.(gb1,Ref(gb2)))
sum(normal_form.(gb2,Ref(gb1)))

normal_form(bg1,gens)

lm.(gb1)
s_polys = [ u[i,2] * rinj(k,j) - u[i,2] * u[k,2] * u[j,3] for i in 2:n, k in 2:n, j in 2:n if i != k && k!=j]

    normal_form.(s_polys,Ref(gens))

lm.(normal_form.(s_polys,Ref(gens)))


_, dct = normal_form_with_rep(bg1,gens)

rinjs
rwels_min
temp = bg1 - rinjs[1]#+rwels_min[1]-rinjs[2]
rinjs[1]
_, dct = normal_form_with_rep(temp,gens)
dct
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

