
#Example: The Quantum Permutation Group on 4 Elements

#= To get this to run

using Pkg
Pkg.add("Oscar#master")
Pkg.add("https://github.com/dmg-lab/QuantumAutomorphismGroups.jl.git")


=#


using Revise	
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
injs = [u[i,j]*u[k,j] for i in 2:n, j in 2:n, k in 2:n if i != k]

rwels = reshape([rwel(k,j,n) for j in 2:n, k in 2:n if j!=k],(n-1)*(n-2))
rinjs = reshape([rinj(k,j,n) for j in 2:n, k in 2:n if j!=k],(n-1)*(n-2))

bg1s = [bg1(i,j,k,u=u) for i in 2:n, j in 2:n, k in 2:n if (k!=i && k!=j)]
bg2s = [bg2(i,j,k,u=u) for i in 2:n, j in 2:n, k in 2:n if (j!=i && k!=j)]
bg3s = [bg3(k,j,u=u) for k in 2:n, j in 2:n if j!=k]
bg4s = [bg4(k,j,u=u) for k in 2:n, j in 2:n if k!=j]
bg5s = [bg5(k,j,u=u) for k in 2:n, j in 2:n if k!=j]
bg6s = [bg6(k,j,u=u) for k in 2:n, j in 2:n if k!=j]
bg7s = [bg7(j,k,h,u=u) for j in 2:n, k in 3:n, h in 2:n if k!=h]
bg8s = [bg8(k,j,h,u=u) for k in 2:n, j in 2:n, h in 2:n if (k!=j && h!=3)]
bg9s = [bg9(k,j,u=u) for k in 3:n, j in 2:n if 3!=j]    # 16
bg10s = [bg10(j,k,u=u) for j in 3:n, k in 2:n if 3!=k]  # 16
bg11s = [bg11(k,i,j,u=u) for k in 3:n, i in 2:n, j in 2:n if i!=j]
bg12s = [bg12(k,j,h,u=u) for k in 2:n, j in 2:n, h in 2:n if (3!=j && k!=j)]
bg13s = [bg13(i,j,k,u=u) for i in 3:n, j in 2:n, k in 2:n if i!=k]
bg14s = [bg14(k,i,j,u=u) for k in 2:n, i in 2:n, j in 2:n if (i!=j && k!=i)]

bg1s_rG0 = [normal_form_noRedTail(f,G0) for f in bg1s]
bg2s_rG0 = [normal_form_noRedTail(f,G0) for f in bg2s]
bg3s_rG0 = [normal_form_noRedTail(f,G0) for f in bg3s]
bg4s_rG0 = [normal_form_noRedTail(f,G0) for f in bg4s]
bg5s_rG0 = [normal_form_noRedTail(f,G0) for f in bg5s]
bg6s_rG0 = [normal_form_noRedTail(f,G0) for f in bg6s]
bg7s_rG0 = [normal_form_noRedTail(f,G0) for f in bg7s]
bg8s_rG0 = [normal_form_noRedTail(f,G0) for f in bg8s]
bg9s_rG0 = [normal_form_noRedTail(f,G0) for f in bg9s]
bg10s_rG0 = [normal_form_noRedTail(f,G0) for f in bg10s]
bg11s_rG0 = [normal_form_noRedTail(f,G0) for f in bg11s]
bg12s_rG0 = [normal_form_noRedTail(f,G0) for f in bg12s]
bg13s_rG0 = [normal_form_noRedTail(f,G0) for f in bg13s]
bg14s_rG0 = [normal_form_noRedTail(f,G0) for f in bg14s]

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

