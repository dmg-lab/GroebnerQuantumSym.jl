
using Oscar.AbstractAlgebra.Generic

export reps_vector_to_poly, reduction_string



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



function reduction_string(
  g1::NamedGenerators,
  ele::Generic.FreeAssociativeAlgebraElem;
  extra::Union{NamedGenerators, Nothing} = nothing,
  to_file::String = "")

  G1, names = generators(g1), g1.names
  


  if !isnothing(extra)
    extra_gens, extra_names = generators(extra), extra.names
    r,v = normal_form_with_rep(ele, extra_gens)

    if iszero(r) 
      @warn "Reduced to zero using only the extra generators"
      return reps_vector_to_poly(v, extra_names)
    end
    y = reps_vector_to_poly(v, extra_names)
    println(y)
  else
    y = ""
  end

  r, v = normal_form_with_rep(ele,G1)
  r != 0 && @warn "This does not reduce to zero, using the normal_form algorithm"
  x = reps_vector_to_poly(v,names)
  
  if y == ""
    x = x
  else
    x = x * " + " * y
  end
  
  if to_file != ""
    open(to_file, "w") do io
        println(io, x)
    end
  end
  return x
end



#= Reduction of rwel23
n = 8
G1 = g1_named(n)
u = magic_unitary(n)
E1 = extra_relations(u)
QuantumGB.add!(E1, G1)

rwel23 = rwel(2,3; u=u)
reduction_string(E1, rwel23)
reduction_string(E1, rwel23; to_file="../data/reduction_strings/n_8_rwel23.txt")
=#

#= Trying todo the same with bg(9,v,h) for v = 5, h = 6
n = 8
G1 = g1_named(n)
u = magic_unitary(n)

#Thing to reduce
v = 5
h = 6
bege9 = bg(9,v,h; u=u)


#Define helper:
es = [ sum([rinj(v,k; u=u) * u[i,h]  for i in 4:n]) for k = 2:n if k != v];
es_names = indexed_name("es", [parse(Int, "$k") for k = 2:n if k != v])
es_ids = [Symbol("es$k") for k in 2:n if k != v]
E1 = named_generators(es, es_names, es_ids)




add!(E1, G1; check=false)


reduction_string(E1, bege9)
reduction_string(G1, bege9; to_file="../data/reduction_strings/n_8_bege9.txt")
reduction_string(G1, bege9)
=#


#=Those are useless (lemma 12)
iter = [(a,k,j,x,y) for a = 2:n for k = 2:n for j = 2:n for x = 3:n for y = 3:n if a != k && k!= j];
srinj = [u[a,2] * rinj(k,j,x,;u=u) for (a,k,j,x,y) in iter];
srinj_names = indexed_name("srinj", [parse(Int, "$a$k$j$x$y") for (a,k,j,x,y) in iter]);
srinj_ids = [Symbol("srinj$a$k$j$x$y") for (a,k,j,x,y) in iter];
add!(E1, srinj, srinj_names, srinj_ids)



n = 8
G0 = g0_named(n)

u = magic_unitary(n)

x  = bg(9,3,4,u=u)+rinj(3,2,u=u)*sum(u[4:n,4]) +sum([rinj(3,i,u=u)*sum(u[2:n,4]) for i in (4,n)]) -u[3,2]*sum([u[i,j] for i in (3:n) for j in (3:n) if j != 4])*col_sum(4,u) +wel(3,2,3,u=u)*sum(u[2:n,4]) -u[3,2]*sum([inj(i,4,j,u=u) for i in (2:n) for j in (2:n) if i != j]) +sum([wel(3,2,i,u=u)*u[j,4] for i in (4:n) for j in (2:n)]) -u[3,2]*sum([ip(i,4,u=u) for i in (3:n)]) -u[3,2]*sum([rwel(i,4,u=u) for i in (5:n)]) +sum([u[3,i]*row_sum(j,u)*u[k,4] for i in (3:n) for j in (2:n) for k in (3:n) if j!=3]) +u[3,3]*sum([row_sum(i,u) for i in (4:n)])*u[2,4] -u[3,3]*sum(u[3:n,2])*col_sum(4,u) +rinj(3,5,u=u)*sum([u[i,4] for i in (2:n)]) +rinj(3,6,u=u)*sum([u[i,4] for i in (2:n)]) +rinj(3,7,u=u)*sum([u[i,4] for i in (2:n)]) -u[3,3]*sum([rwel(i,4,u=u) for i in (2:n) if i!=3 && i!=4]) +wel(3,3,2,u=u)*sum([u[i,4] for i in (2:n)]) -inj(3,3,2,u=u)*sum([u[i,4] for i in (3:n)]) -sum([inj(3,3,j,u=u)*u[i,4] for i in (2:n) for j in (4:n)]) -u[3,3]*sum([inj(j,4,i,u=u) for i in (2:n) for j in (2:n) if i != j]) +sum([wel(3,3,j,u=u)*u[i,4] for i in (2:n) for j in (4:n)]) -u[3,3]*sum([ip(i,4,u=u) for i in (3:n)]) -u[3,3]*sum([u[i,j] for i in (3:n) for j in (5:n)])*col_sum(4,u)

reduction_string(G0, x)
r,v = normal_form_with_rep(x, G0.gs)
rs4 = G0[:rs4]
r1,v1 = normal_form_with_rep(x,[rs4]);

ok, ml, mr = Generic.word_divides_leftmost(x.exps[1], rs4.exps[1])
v
x

=#
