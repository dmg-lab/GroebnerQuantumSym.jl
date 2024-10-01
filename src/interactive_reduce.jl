
using Oscar.AbstractAlgebra.Generic

export reps_vector_to_poly



"""
```jldoctest
using Oscar
R, (x, y) = free_associative_algebra(QQ, ["x", "y"]);
f = 2x^2 + 3x*y*x^2 + 4y^2;
g = [x^2, y^2, x*y];

zr,  dct = normal_form_with_rep(f, g);
reps_vector_to_poly(dct)
# output

Dict{AbstractAlgebra.Generic.FreeAssociativeAlgebraElem, AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{QQFi
eldElem}} with 2 entries:
  x^2 => 2*x^2
  y^2 => 4*y^2
```

"""
function reps_vector_to_poly(v::Vector)
  length(v) == 0 && return 0
  A = parent(v[1][3])
  f = zero(A)
  for rep in v
    f += rep[1] * rep[2] * rep[3] * rep[4]
  end
  return f
end
"""
```jldoctest
using Oscar
R, (x, y) = free_associative_algebra(QQ, ["x", "y"]);
f = 2x^2 + 3x*y*x^2 + 4y^2;
g = [x^2, y^2, x*y];
names = Dict(x^2 => "u", y^2 => "v", x*y => "w")

zr,  dct = normal_form_with_rep(f, g);
reps_vector_to_poly(dct, names)
# output

Dict{AbstractAlgebra.Generic.FreeAssociativeAlgebraElem, AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{QQFi
eldElem}} with 2 entries:
  x^2 => 2*x^2
  y^2 => 4*y^2
```

"""
function reps_vector_to_poly(v::Vector, names::Dict{<:Generic.FreeAssociativeAlgebraElem,String})
  v3 = [x[3] for x in v]
  @assert all([x in keys(names) for x in v3])
  length(v) == 0 && return 0
  f = ""
  A = parent(v[1][3])
  for rep in v
    if rep[1] != one(A)
      f *= string(rep[1])
      f *= "*"
    end
    if rep[2] != one(A) 
      f *= string(rep[2]) 
      f *= "*"
    end
    f *= names[rep[3]] 
    if rep[4] != one(A) 
      f *= "*"
      f *= string(rep[4])
    end
    f *= " + "
  end
return f[1:end-2]
end

#= The Example

n = 4
G0, names = g0_extended(n,names=true)
u = magic_unitary(G0)
bg1 = rwel(2,3; u = u)
r, v = normal_form_with_rep(bg1,G0)
r != 0 && error("r != 0")
x = reps_vector_to_poly(v,names)
length(x)
print(x)

n = 8
G0, names = g0(n,names=true)
u = magic_unitary(G0)
bg1 = bg(1,2,3,4; u = u)
r, v = normal_form_with_rep(bg1,G0)
r != 0 && error("r != 0")
x = reps_vector_to_poly(v,names)
r
groebner_basis(G0)
length(x)
print(x)



n = 4
G0, names = g0(n,names=true);
g1 = groebner_basis(G0)
x = leading_monomial.(b1)
showall(x)
for x in gb1
open("../examples/gb4.txt", "a") do io
    println(io, x)
end
end



u = magic_unitary(G0);



using Oscar
n = 6
G0, names = g0(n,names=true);
G1,names = g1(n,true)
G0, names = g0_extended(n,names=true);


u = magic_unitary(G0)
x1 = sum([u[j,1] - u[k,2] - u[1,j] + u[2,k] for k in 2:n for  j in 2:n if j != k])

normal_form(x1,G0)

r, v = normal_form_with_rep(x1,G0)
r != 0 && error("r != 0")
x = reps_vector_to_poly(v,names)
print(x)

rwel(2,3; u = u) in G0
groebner_basis(G0)

rinj(2,3; u = u)
gb1 = groebner_basis(G0,6)
u = magic_unitary(G0)

filter(x -> x != 0, red_bgs)


[x in gb1 for x in bg1s]


r, v = normal_form_with_rep(bg1s[1],G0)
v
r != 0 && error("r != 0")
x = reps_vector_to_poly(v,names)
r
length(x)
print(x)

n = 6
G1,names = g1(n,true)
u = magic_unitary(G1)

o1=bg(2,3,6,3; u = u)
o2=bg(8,3,6,2; u = u)

o2*u[6,4]*u[3,3]
u[2,3]*u[4,6]*o1

ov1 = o2*u[6,4]*u[3,3] -u[2,3]*u[4,6]*o1

r, v = normal_form_with_rep(ov1,G1)
r != 0 && error("r != 0")
x = reps_vector_to_poly(v,names)


rinj46_1 = rinj(4,2; u = u)u[3,6] -u[4,2]u[2,3]u[3,6]
rinj46_2 = u[4,2]rwel(3,6; u = u) -u[4,2]u[2,3]u[3,6]
bege1 = bg(9,4,6; u = u)

rinj46_1-rinj46_2 == bege1 #should be true

r, v = normal_form_with_rep(rinj46,G1)

r, v = normal_form_with_rep(bege1,G1)
r != 0 && error("r != 0")
x = reps_vector_to_poly(v,names)
open("./bege1.txt", "w") do io
    println(io, x)
end
print(x)


=#
