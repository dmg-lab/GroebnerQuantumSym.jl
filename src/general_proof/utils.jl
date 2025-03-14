
export Lexicon, lexicon, Term, create_preimage_in


function reverse_dict(d::Dict)
  Dict(v=>k for (k,v) in d)
end

function unique_eles(v::Vector{Vector{T}}) where T
  s = Set{T}()
  for x in v
    for y in x
      push!(s, y)
    end
  end
  return s
end

mutable struct Lexicon
  predicates::Vector{Vector{Int}}
  labels::Dict{String,Int}
  coefficient_ring::Ring
  function Lexicon(v::Vector{Vector{Int}}, dct::Dict{String,Int}, coeff_ring::Ring=QQ[:n][1])
    @assert all(x->length(x) == length(v[1]), v) "Not all elements in v have the same length"
    words = unique_eles(v)
    @assert all(x->haskey(reverse_dict(dct), x), words) "Not all elements in v are in the dictionary"
    new(v, dct, coeff_ring)
  end
end

function Base.show(io::IO, l::Lexicon)
  println(io, "Lexicon with $(length(l.predicates)) predicates")
end


"""
  lexicon(v::Vector{Vector{Int}})

Basic constructor for a lexicon. 

Example:

```julia
L1 = lexicon([[1],[2],[3],[-4]])
L2 = lexicon([[1],[2],[3],["g"]])

```
"""
function lexicon(v::Vector{Vector{Int}}, dct::Dict{String,Int}, coeff_ring::Ring=QQ[:n][1])
  @assert all(x->length(x) == length(v[1]), v) "Not all elements in v have the same length"
  @assert all(all(x->haskey(reverse_dict(dct), x), y) for y in v) "Not all elements in v are in the dictionary"
  Lexicon(v, dct, coeff_ring)
end

function lexicon(v::Vector{Vector{Int}})
  @assert all(x->length(x) == length(v[1]), v)
  dct = Dict{String,Int}()
  # set of all integers in the lexicon
  words = unique_eles(v)
  for x in words
    if x > 0
      dct[string(x)] = x
    else
      dct["≥$(abs(x))"] = x
    end
  end
  lexicon(v, dct)
end


function lexicon(v::Vector{Vector{String}})
  @assert all(x->length(x) == length(v[1]), v)
  wrds = unique_eles(v)
  new_v = Vector{Int}[]
  dct = Dict{String,Int}()
  for (i, x) in enumerate(wrds)
    dct[x] = i
  end
  for x in v
    push!(new_v, [dct[y] for y in x])
  end
  lexicon(new_v, dct)
end

function lexicon(v::Vector{Vector})
    v = Vector{String}[string.(x) for x in v]
    lexicon(v)
end

function find(l::Lexicon, pred::Vector{Int})
  findfirst(x->x==pred, l.predicates)
end

"""
  Base.:*(l1::Lexicon, l2::Lexicon)

This function computes the cartesian product of two lexicons.

Example:
```julia
L1 = lexicon([[1],[2],[3],["g"]])
L2 = lexicon([[1],[2],["g"],[3]])

L = L1 * L2
```
"""
function Base.:*(l1::Lexicon, l2::Lexicon)
  new_rep = Dict{String,Int}()
  @assert l1.coefficient_ring == l2.coefficient_ring "Coefficient rings do not match"

  rv1 = reverse_dict(l1.labels)
  rv2 = reverse_dict(l2.labels)
  i = 1
  for (k,v) in l1.labels
    new_rep[k] = i
    i += 1
  end
  for (k,v) in l2.labels
    if !haskey(new_rep, k)
      new_rep[k] = i
      i += 1
    end
  end

  new_v = Vector{Int}[]
  for x in l1.predicates
    for y in l2.predicates
      reps_x = [new_rep[rv1[z]] for z in x]
      reps_y = [new_rep[rv2[z]] for z in y]
      push!(new_v, vcat(reps_x, reps_y))
    end
  end

  lexicon(new_v, new_rep, l1.coefficient_ring)
end

mutable struct Term
  L::Lexicon
  pairs::Vector{Tuple{RingElem,Int}}
  function Term(L::Lexicon, pairs::Vector{Tuple{RingElem,Int}})
    @assert (all(x->x[2] <= length(L.predicates), pairs)) "Index out of bounds"
    @assert (all(x->parent(x[1]) == L.coefficient_ring, pairs)) "Coefficents not in the coefficient ring"
    new(L, pairs)
  end
  Term(L::Lexicon) = new(L, [])
end

function Base.push!(t::Term, x::Tuple{RingElem,Int})
  @assert x[2] <= length(t.L.predicates) "Index out of bounds"
  @assert parent(x[1]) == t.L.coefficient_ring "Coefficient not in the coefficient ring"
  push!(t.pairs, x)
end

Base.iszero(t::Term) = isempty(t.pairs)

function Base.show(io::IO, t::Term)
  if iszero(t)
    println("Zero term")
    return
  end

  rev_dct = reverse_dict(t.L.labels)
  ret = map(x->(x[1], [rev_dct[y] for y in t.L.predicates[x[2]]]), t.pairs) 
  println("Term over lexicon with $(length(t.L.predicates)) predicates:")
  show(io, "text/plain", ret)
end


function create_preimage_in(L::Lexicon, v::Vector{Vector{Any}}, coeff::RingElem; filter::Vector = [])
  v = Vector{String}[string.(x) for x in v]
  s = Iterators.product(v...) # ("str", "str", "str")...
  labels = L.labels
  P = Vector{Int}[]
  for (i, j) in filter
    if j > i
      i, j = j, i
    end
    s = Iterators.filter(x -> x[i] == "g"  || x[j] == "g" || x[i] != x[j], s)
  end 

  for x in s
    x = collect(x) 
    for (i, j) in filter
      if x[j] == "g" && x[i] == "g"
        x[i] = "g⩓≠i$(j)"
      end
    end
    push!(P, [labels[y] for y in x])
  end

  vec = Term(L)
  R = L.coefficient_ring
  for x in P  
    idx = find(L, x)
    @assert !isnothing(idx) "$x not found in Lexicon"
    push!(vec, (R(coeff), idx))
  end
  return vec
end

function create_preimage_in(L::Lexicon, v::Vector{Vector{Any}}, coeff::Number; filter::Vector = [])
  R = L.coefficient_ring
  coeff = R(coeff)
  return create_preimage_in(L, v, coeff, filter=filter)
end

function create_preimage_in(L::Lexicon, v::Vector{Vector}, coeff::RingElem; filter::Vector = [])
  v = Vector{Any}[string.(x) for x in v]
  return create_preimage_in(L, v, coeff, filter=filter)
end

function create_preimage_in(L::Lexicon, v::Vector{Vector}, coeff::Number; filter::Vector = [])
  R = L.coefficient_ring
  coeff = R(coeff)
  return create_preimage_in(L, v, coeff, filter=filter)
end


function Base.:+(t1::Term, t2::Term)
  @assert t1.L == t2.L "Lexicons do not match"
    new_pairs = copy(t1.pairs)
  for x in t2.pairs
    i = findfirst(y->y[2] == x[2], new_pairs)
    if isnothing(i)
      push!(new_pairs, x)
    else
      if new_pairs[i][1] + x[1] == 0
        deleteat!(new_pairs, i)
      else
        new_pairs[i] = (new_pairs[i][1] + x[1], new_pairs[i][2])
      end
    end
  end
  return Term(t1.L, new_pairs)
end

