
#Make it more readable should probably be changed tho.
Base.getindex(expr::Expr, i::Int) = expr.args[i]
Base.setindex!(expr::Expr, value, i::Int) = (expr.args[i] = value; expr)


#= Section IndexSet =#

mutable struct IndexSet
  var::Symbol
  start_bound::Union{Int, Symbol}
  end_bound::Union{Int, Symbol}
end

function Base.show(io::IO, idx::IndexSet) 
  println(io, "$(idx.var) ∈ {$(idx.start_bound):$(idx.end_bound)}")
end

#= Section Normal Form =#

mutable struct NF 
  degree::Int
end

#= Example
  code = "u[a,2]*u[c,d]"
  expr = Meta.parse(code)
  _get_degree(expr)
=#
function _get_degree(expr::Expr)
  if expr.head === :ref
    @assert expr.args[1] === :u "Only :u is allowed"
    return 1, [expr.args[2],expr.args[3]]
  elseif  expr.head === :call && expr.args[1] === :*
    d1, s1 = _get_degree(expr.args[2])
    d2, s2 = _get_degree(expr.args[3])
    return d1+d2, vcat(s1,s2)
  else
    @assert false "Invalid expression: $(expr)"
  end
end

function _is_prod_of_refs(expr::Expr)
  if expr.head === :ref
    return true
  elseif expr.head === :call && expr.args[1] === :*
  return all(x-> _is_prod_of_refs(x), expr.args[2:end])
  else
    return false
  end
end



#= Section IntermediateRepr =#

mutable struct IntermediateRepr
  idx_sets::Vector{IndexSet}
  cond_set::Vector{Expr}
  formula::Union{Nothing,Expr,NF}
  ItermediateRep() = new(IndexSet[], Expr[], nothing)
end

function normalize!(interrep::IntermediateRepr)
  if interrep.formula isa NF
    println("Already normalized")
    return interrep
  elseif interrep.formula isa Expr
    if _is_prod_of_refs(interrep.formula)
      degree, idxs = _get_degree(interrep.formula)
      #Add constants
      for (i,ids) in enumerate(idxs)
        if ids isa Int
          ids_new = Symbol("i$i")
          push!(interrep.idx_sets, IndexSet(ids_new, ids, ids))
          ids = ids_new
        end
      end


      # Resolve the idxs_sets
      for idx_set in interrep.idx_sets
        i = findfirst(x->x==idx_set.var, idxs) 
        isnothing(i) && continue
        idx_set.var = Symbol("i$i")
      end
      #Resolve the condition set THIS IS YANKY AF
      new_cond_set = Expr[]
      sizehint!(new_cond_set, length(interrep.cond_set))
      for cond in interrep.cond_set
        cond_str = string(cond)
        for (i,idx) in enumerate(idxs)
          cond_str = replace(cond_str, string(idx) => "i$i")
        end
        push!(new_cond_set, Meta.parse(cond_str))
      end
      interrep.cond_set = new_cond_set
      interrep.formula = NF(degree)
    end
  else
    @assert false "Empty formula"
  end
  return interrep
end



#= 
code = "[u[i,j] for i in 2:n for j in 5:n if i != j]"
expr = Meta.parse(code)

Ir = resolve_condition!(expr)
normalize!(Ir)

code1 = "[u[i,2]*u[k,h] for i in 2:n for k in 1:n for h in 1:n if i != h && k != h]"
Ir1 = resolve_condition!(Meta.parse(code1))
normalize!(Ir1)

code2 = "[u[1,j] for j in 2:n]"
Ir2 = resolve_condition!(Meta.parse(code2))
normalize!(Ir2)

tst = Meta.parse("u[1,2]*u[3,4]*u[5,6]*sum([u[i,j] for i in 2:n for j in 5:n if i != j])")
_is_prod_of_refs(tst)


=#

