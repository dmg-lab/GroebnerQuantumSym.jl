
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

n = 8
G0, names = g0_extended(n,names=true)
u = magic_unitary(G0)
bg1 = rwel(2,3; u = u)
r, v = normal_form_with_rep(bg1,G0)
r != 0 && error("r != 0")
x = reps_vector_to_poly(v,names)
length(x)
print(x)

G0, names = g0(n,names=true);
u = magic_unitary(G0);
bg1 = rwel(2,3; u = u);
r, v = normal_form_with_rep(bg1,G0);
r != 0 && error("r != 0");
x1 = reps_vector_to_poly(v,names);
length(x1)
print(x1)



n = 4
G0, names = g0(n,names=true);
gb1 = groebner_basis(G0)
for x in gb1
open("../examples/gb4.txt", "a") do io
    println(io, x)
end
end



u = magic_unitary(G0);



using Oscar
n = 4
G0, names = g0(n,names=true);
gb1 = groebner_basis(G0)
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

# For multiple n,
using Oscar
for n in 11:25
  #Check if already in the file
  G0, names = g0(n,names=true)
  u = magic_unitary(G0)
  bg1 = rwel(2,3; u = u)
  r, v = normal_form_with_rep(bg1,G0)
  r != 0 && error("r != 0")
  x = reps_vector_to_poly(v,names)
  x = "n = $(n): " * x
  #print to file, or append to file
  open("../examples/rwel_23_reductions.txt", "a") do io
    println(io, x)
  end
end

print(x)
names

=#
