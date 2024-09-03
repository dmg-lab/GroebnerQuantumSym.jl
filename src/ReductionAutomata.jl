#An Auntomata for interactive reduction of a given polynomial
import Base.delete!


export ReductionStep,
  NamedGenerators,
  div_by,
  reduction_polynomial,
  new_f,
  monomial_at_position,
  named_generators,
  delete!,
  add!,
  reduction_automata,
  generators
  polynomial,
  current_state,


"""
  monomial_at_position(f::Generic.FreeAssociativeAlgebraElem{T}, i::Int)

Given a polynomial `f` and an integer `i`, return the monomial at position `i`.

```julia
using Oscar
R, (x, y) = free_associative_algebra(QQ, ["x", "y"]);
f = 2x^2*y*x + 3x*y + 4y^2;

monomial_at_position(f, 2)

```
"""
function monomial_at_position(f::Generic.FreeAssociativeAlgebraElem{T}, i::Int) where {T<:FieldElement}
    @assert length(f) >= i "f does not have enough monomials"
    A = parent(f)
    c = f.coeffs[i]
    mr = Int[]
    return Generic.mul_term(c, f.exps[i], one(A), mr)
end

function Base.getindex(f::Generic.FreeAssociativeAlgebraElem{T}, i::Int) where {T<:FieldElement}
    return monomial_at_position(f, i)
end



"""
  ReductionStep

A struct to store the information of a reduction step.

```julia
using Oscar
R, (x, y) = free_associative_algebra(QQ, ["x", "y"]);
f = 2x^2*y*x + 3x*y + 4y^2;
g = x

div_by(f, g)

# output

3*x*y + 2*y^2
```

"""
mutable struct ReductionStep
    old_f::Generic.FreeAssociativeAlgebraElem{T} where {T<:FieldElement}
    c::FieldElement
    ml::Vector{Int}
    g::Generic.FreeAssociativeAlgebraElem{T} where {T<:FieldElement}
    mr::Vector{Int}
    function ReductionStep(old_f::Generic.FreeAssociativeAlgebraElem{T},
      c::T,
      ml::Vector{Int},
      g::Generic.FreeAssociativeAlgebraElem{T},
      mr::Vector{Int}) where {T<:FieldElement}
      @assert parent(old_f) == parent(g)
        new(old_f, c, ml, g, mr)
    end
end

function div_by(f::Generic.FreeAssociativeAlgebraElem{T}, g::Generic.FreeAssociativeAlgebraElem{T}; right::Bool=false, i::Int=1, j::Int=1) where {T<:FieldElement}
    @assert length(f.exps) >= i "f does not have enough monomials"
    @assert length(g.exps) >= j "g does not have enough monomials"

    if right
        ok, ml, mr = Generic.word_divides_rightmost(f.exps[i], g.exps[j])
    else
        ok, ml, mr = Generic.word_divides_leftmost(f.exps[i], g.exps[j])
    end
    @assert ok # we assume that g divides f
    if ok
        qi = divexact(f.coeffs[i], g.coeffs[j])
        return ReductionStep(f,  qi, ml, g, mr)
    end
    error("monomial at position $i does not divide")
end

"""
  reduction_polynomial(R::ReductionStep)

Given a reduction step `R`, return the polynomial used to reduce `R.old_f`.

```julia
using Oscar
R, (x, y) = free_associative_algebra(QQ, ["x", "y"]);
f = 2x^2*y*x + 3x*y + 4y^2;
g = x

left = div_by(f, g, right=false)
right = div_by(f, g, right=true)

reduction_polynomial(left)
reduction_polynomial(right)
```

"""
function reduction_polynomial(R::ReductionStep)
  return Generic.mul_term(R.c, R.ml, R.g, R.mr)
end

"""
  new_f(R::ReductionStep)

Given a reduction step `R`, return the new polynomial after the reduction.

```julia

R, (x, y) = free_associative_algebra(QQ, ["x", "y"]);
f = 2x^2*y*x + 3x*y + 4y^2;
g = x

Rs = div_by(f, g)
new_f(Rs)

```
"""
function new_f(R::ReductionStep)
  return Generic._sub_rest(R.old_f, reduction_polynomial(R), 1)
end

function Base.show(io::IO, R::ReductionStep)
    otp = ""
    otp *= "Reduction: f --- $(reduction_polynomial(R)) ---> $(new_f(R))\n"
    print(io, otp)
end

"""
  NamedGenerators

A struct to store the names of the generators of an Ideal in a Free Associative Algebra.

```julia
using Oscar
R, (x, y) = free_associative_algebra(QQ, ["x", "y"]);
gs = [x^2, y^2, x*y];
names = Dict(x^2 => "g₁", y^2 => "g₂", x*y => "g₃")
ids = Dict(:g1 => x^2, :g2 => y^2, :g3 => x*y)

ng = named_generators(gs, names, ids)

delete!(ng, y^2)
delete!(ng, :g1)

ng
ng[ng[:g1]]
```
"""
mutable struct NamedGenerators
    gs::Vector{Generic.FreeAssociativeAlgebraElem{T}} where {T<:FieldElement}
    names::Dict{<:Generic.FreeAssociativeAlgebraElem, String}
    identifiers::Dict{Symbol, Int}
end

function named_generators(
  gs::Vector{Generic.FreeAssociativeAlgebraElem{T}},
  names::Dict{Generic.FreeAssociativeAlgebraElem{T}, String},
  ids::Dict{Symbol, Int}) where {T<:FieldElement}

  n = length(gs)
  @assert n == length(names) == length(values(ids))
  @assert all([x in keys(names) for x in gs])
  @assert all([x <= n  for x in values(ids)])

  return NamedGenerators(gs, names, ids)
end

function named_generators(
  gs::Vector{Generic.FreeAssociativeAlgebraElem{T}},
  names::Dict{Generic.FreeAssociativeAlgebraElem{T}, String},
  ids::Dict{Symbol, Generic.FreeAssociativeAlgebraElem{T}}) where {T<:FieldElement}
  
  @assert length(gs) == length(values(ids))
  @assert all([x in gs for x in values(ids)])
  
  identifiers = Dict{Symbol, Int}()
  for (k, v) in ids
    identifiers[k] = findfirst(x -> x == v, gs)
  end
  return named_generators(copy(gs), copy(names), copy(identifiers))
end

function _get_ids(ng::NamedGenerators)
  return Dict(value => key for (key, value) in ng.identifiers)
end


Base.getindex(ng::NamedGenerators, i::Int) = ng.gs[i]
Base.setindex!(ng::NamedGenerators, g::Generic.FreeAssociativeAlgebraElem, i::Int) = ng.gs[i] = g
Base.length(ng::NamedGenerators) = length(ng.gs)
Base.getindex(ng::NamedGenerators, s::Symbol) = ng.gs[ng.identifiers[s]]
Base.setindex!(ng::NamedGenerators, g::Generic.FreeAssociativeAlgebraElem, s::Symbol) = ng.gs[ng.identifiers[s]] = g
Base.getindex(ng::NamedGenerators, f::Generic.FreeAssociativeAlgebraElem) = ng.names[f]
Base.setindex!(ng::NamedGenerators, name::String, f::Generic.FreeAssociativeAlgebraElem) = ng.names[f] = name



function delete!(ng::NamedGenerators, f::Generic.FreeAssociativeAlgebraElem{T}) where {T<:FieldElement}
  index = findfirst(x -> x == f, ng.gs)

  delete!(ng.names, f)
  rev_ids = QuantumGB._get_ids(ng)
  println(rev_ids)
  delete!(ng.identifiers, rev_ids[index])
  deleteat!(ng.gs, index)
  #Reduce the identifiers after the position of the deleted element
  for (k, v) in ng.identifiers
    if v > index
      ng.identifiers[k] -= 1
    end
  end
  return ng 
end

delete!(ng::NamedGenerators, s::Symbol) = delete!(ng, ng[s])

Base.show(io::IO, ng::NamedGenerators) = print(io, "Generating set with $(length(ng)) elements")

function add!(ng::NamedGenerators, f::Generic.FreeAssociativeAlgebraElem{T}, name::String, id::Symbol) where {T<:FieldElement}
  push!(ng.gs, f)
  ng.names[f] = name
  ng.identifiers[id] = length(ng.gs)
  return ng
end

function add!(ng::NamedGenerators,
  fs::Vector{Generic.FreeAssociativeAlgebraElem{T}},
  names::Vector{String},
  ids::Vector{Symbol}) where {T<:FieldElement}
  @assert length(fs) == length(names) == length(ids)
  for i in 1:length(fs)
    add!(ng, fs[i], names[i], ids[i])
  end
end

function named_generators(
  fs::Vector{Generic.FreeAssociativeAlgebraElem{T}},
  names::Vector{String},
  ids::Vector{Symbol}) where {T<:FieldElement}
  ng = NamedGenerators(Vector{Generic.FreeAssociativeAlgebraElem{T}}(), Dict{Generic.FreeAssociativeAlgebraElem{T}, String}(), Dict{Symbol, Int}())
  add!(ng, fs, names, ids)
  return ng
end

generators(ng::NamedGenerators) = ng.gs

mutable struct ReductionAutomata
    generator_set::NamedGenerators
    f::Generic.FreeAssociativeAlgebraElem{T} where {T<:FieldElement}
    steps::Vector{ReductionStep} 
    A::FreeAssociativeAlgebra
    magic_unitary::Matrix{Generic.FreeAssociativeAlgebraElem{T}} where {T<:FieldElement}
    current::Generic.FreeAssociativeAlgebraElem{T} where {T<:FieldElement}
end

function check_integrity(Rs::ReductionAutomata)
  @assert all([parent(x.old_f) == parent(Rs.f) for x in Rs.steps])
  # More to come
  return nothing
end

"""
  reduction_automata(generator_set::NamedGenerators, f::Generic.FreeAssociativeAlgebraElem{T})

Given a named generator set and a polynomial `f`, return a `ReductionAutomata` object.

```julia
using Oscar
G0 = g0_named()
u = magic_unitary(G0)
f = rwel(2,3; u=u)

Rs = reduction_automata(G0, f)
```
"""
function reduction_automata(
  generator_set::NamedGenerators,
  f::Generic.FreeAssociativeAlgebraElem{T}) where {T<:FieldElement}
  Rs = ReductionAutomata(generator_set, f, ReductionStep[], parent(f),magic_unitary(generator_set), f)
  check_integrity(Rs)
  return Rs
end

Base.show(io::IO, Rs::ReductionAutomata) = print(io, "Reduction Automata with $(length(Rs.steps)) steps over $(length(Rs.generator_set)) generators")

function finish!(Rs::ReductionAutomata)
  error("Not implemented")
end

polynomial(Rs::ReductionAutomata) = Rs.f
current_state(Rs::ReductionAutomata) = Rs.current
