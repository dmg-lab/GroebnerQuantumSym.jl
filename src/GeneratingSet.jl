
export name,add,delete

mutable struct GeneratingSet
  generators::Vector{Generic.FreeAssAlgElem}
  parent::Generic.FreeAssAlgebra
  u::Array{Generic.FreeAssAlgElem,2}
  names::Dict{Generic.FreeAssAlgElem,String}
  function GeneratingSet(n::Int)
    relat , relat_by_type, u, A = getQuantumPermutationGroup(n);
    ans = new()
    ans.generators = relat
    ans.parent = A
    ans.u = u
    return ans
  end
  function GeneratingSet(named_gens::Dict{Generic.FreeAssAlgElem{T},String}) where T
    @assert length(named_gens) > 0
    ans = new()
    ans.generators = collect(keys(named_gens))
    prt = parent(ans.generators[1])
    ans.parent = prt
    ans.names = named_gens
    n = Int(sqrt(length(prt.S)))
    u = Array{Generic.FreeAssAlgElem,2}(undef, n, n)
    for i in 0:(n-1)
      for j in 1:n
        u[i+1,j] = prt[i*n+j]
      end
    end
    ans.u = u
    return ans
  end
end

function name(x::Generic.FreeAssAlgElem,genset::GeneratingSet) 
  @assert x in genset.generators "Element not in generating set"
  @assert x in keys(genset.names) "Element has no name in this generating set"
  genset.names[x]
end

function relation(str::String, genset::GeneratingSet)
  @assert str in values(genset.names) "Relation not a name of a generator"
  for (x,y) in genset.names
    if y == str
      return x
    end
  end
end

Base.getindex(g::GeneratingSet, str::String) = relation(str,g)
Base.getindex(g::GeneratingSet, i::Int) = g.generators[i]
Base.getindex(g::GeneratingSet, x::Generic.FreeAssAlgElem) = name(x,g)
Base.setindex!(g::GeneratingSet, x::Generic.FreeAssAlgElem, str::String) = g.names[x] = str
Base.setindex!(g::GeneratingSet, str::String, x::Generic.FreeAssAlgElem) = g.names[x] = str
Base.show(io::IO, g::GeneratingSet) = print("Generating set with $(length(g.generators)) generators")

Base.length(g::GeneratingSet) = length(g.generators)

function delete(g::GeneratingSet, x::Generic.FreeAssAlgElem)
  @assert x in g.generators "Element not in generating set"
  delete!(g.names,x)
  filter!(y->y != x, g.generators)
end

function add(g::GeneratingSet, x::Generic.FreeAssAlgElem, str::String)
  @assert !(x  in g.generators) "Element already in generating set"
  g.names[x] = str
  push!(g.generators,x)
end

mutable struct ReductionAutomata
  genset::GeneratingSet
  states::Vector{Generic.FreeAssAlgElem}
  finished::Bool
end

function reduction_automata(g::GeneratingSet, f::Generic.FreeAssAlgElem)
  ReductionAutomata(g,[f],false)
end

finish!(r::ReductionAutomata) = r.finished = true
is_finished(r::ReductionAutomata) = r.finished

function next!(R::ReductionAutomata)
  @assert !is_finished(R) "Automata is finished"
  curr_f = R.states[end]
  iszero(curr_f) && finish!(R) && return 

  genset = R.genset
   
  s = length(genset)
    @label first
    i = 1
    @label again 
    if iszero(f)
      println(otp)
      return f 
    end
    ok, ml, mr = Generic.word_divides_leftmost(f.exps[1], g[i].exps[1])
    if !ok && i < s
        i += 1
        @goto again
    end


end





#=
using Oscar
using QuantumAutomorphismGroups
using QuantumGB
x = GeneratingSet(4)
x.generators[1]
dct = Dict(x.generators[1] => "idempotent")
y = GeneratingSet(dct)
y[y["idempotent"]] = "ip_11"
y[x.generators[2]] = "row_sum_1"
add(y,x.generators[3],"col_sum_1")
y.names

name(y.generators[1],y)
y.generators

length(x.parent.S)

=#

