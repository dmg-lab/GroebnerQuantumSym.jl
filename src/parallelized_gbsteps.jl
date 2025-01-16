using Oscar.AbstractAlgebra.Generic

export get_all_divisors,recursive_reduction, normal_form_with_rep_choices


"""
  parallelized_gbsteps(f::Generic.FreeAssociativeAlgebraElem{T}, g::Vector{Generic.FreeAssociativeAlgebraElem{T}}) where {T}

Get all the divisors of the leading monomial f by g.
```jldoctest
using Oscar
R, (x, y) = free_associative_algebra(QQ, ["x", "y"]);
f = 2x^2 + 3x*y*x^2 + 4y^2;
g = [x^2, y^2, x*y];

n = 6
G1 = g1_named(n)
u = magic_unitary(n)
rwel23 = rwel(2,3; u=u)
g = G1.gs
reduction_string(G1,rwel23-G1["rinj23"])
get_all_divisors(rwel23 - G1["rinj23"], g)

```

"""
function get_all_divisors(
  f::Generic.FreeAssociativeAlgebraElem{T},
  g::Vector{Generic.FreeAssociativeAlgebraElem{T}}) where {T}
  
  A = parent(f)
  reps = Tuple{elem_type(base_ring(A)),Generic.FreeAssociativeAlgebraElem{T},Generic.FreeAssociativeAlgebraElem{T},Generic.FreeAssociativeAlgebraElem{T}}[]
  s = length(g)

  iszero(f) && return reps

  for ele in g
    qi = divexact(f.coeffs[1], ele.coeffs[1])
    ok, ml, mr = Generic.word_divides_leftmost(f.exps[1], ele.exps[1])
    if !ok
      continue
    end
    mul_l = one(A)
    mul_r = one(A)

    if length(ml) > 0
      mul_l = prod([A[i] for i in ml])
    end
    if length(mr) > 0
      mul_r = prod([A[i] for i in mr])
    end
    push!(reps, (qi, mul_l, ele, mul_r))
  end
  return reps
end

#=
n = 6
G1 = g1_named(n)
u = magic_unitary(n)
rwel23 = rwel(2,3; u=u)
g = G1.gs
recursive_reduction(rwel23,g)

u = magic_unitary(n)

bg934 = bg(9,3,4,u=u)
recursive_reduction(bg934, G1.gs)
bege8243 = bg(8,2,4,3; u=u)
=#



function recursive_reduction(
  f::Generic.FreeAssociativeAlgebraElem{T},
  g::Vector{Generic.FreeAssociativeAlgebraElem{T}};
  A::Generic.FreeAssociativeAlgebra{T}=parent(f),
  reps=Tuple{elem_type(base_ring(A)),Generic.FreeAssociativeAlgebraElem{T},Generic.FreeAssociativeAlgebraElem{T},Generic.FreeAssociativeAlgebraElem{T}}[],
  base_case=-1,
  open::Int=1,
  logging_length::Int=-1
  ) where {T}

  logging_max = -1

function _recursive_reduction(
  f::Generic.FreeAssociativeAlgebraElem{T},
  g::Vector{Generic.FreeAssociativeAlgebraElem{T}};
  A::Generic.FreeAssociativeAlgebra{T}=parent(f),
  reps=Tuple{elem_type(base_ring(A)),Generic.FreeAssociativeAlgebraElem{T},Generic.FreeAssociativeAlgebraElem{T},Generic.FreeAssociativeAlgebraElem{T}}[],
  base_case=-1,
  open::Int=1,
  logging_length::Int=-1
  ) where {T}
  base_length=-1

  if base_case != -1
    base_length = length(base_case[2])
  end

  iszero(f) && return f, reps

  if base_length != -1 && length(reps) >= base_length - 1
    return f, -1
  end


  choices = get_all_divisors(f, g)



  open += length(choices)
  for (qi, _, ele, _) in choices
        open -= 1
        #Inefficient, but we need to copy the reps
        _, ml, mr = Generic.word_divides_leftmost(f.exps[1], ele.exps[1])
        new_f = Generic._sub_rest(f, Generic.mul_term(qi, ml, ele, mr), 1)
        mul_l = one(A)
        mul_r = one(A)

        if length(ml) > 0
          mul_l = prod([A[i] for i in ml])
        end
        if length(mr) > 0
          mul_r = prod([A[i] for i in mr])
        end
        push!(reps, (qi, mul_l, ele, mul_r))
        #break if computation takes longer then other branch
        r_call = _recursive_reduction(
        new_f, g, A=A, reps=deepcopy(reps), 
         base_case=base_case, open=open)
        #If branch was pruned, continue
        r_call[2] == -1 && continue
        
        if base_length == -1
          base_case = r_call
          base_length = length(base_case[2])
          if logging_max == -1
            println("Found a base case of length $(base_length)")
          elseif base_length < logging_max
            println("Improved the base case from $(logging_max) to $(length(r_call[2]))")
          end
          logging_max = length(r_call[2])
        elseif base_length > length(r_call[2])
        #If we found a new shorter path, set it as the base case
          base_case = r_call
          base_length= length(base_case[2])
          if logging_max == -1
            println("Found a base case of length $(base_length)")
          elseif base_length < logging_max
            println("Improved the base case from $(logging_max) to $(length(r_call[2]))")
          end
            logging_max = length(r_call[2])
        end
  end
  if base_case == -1
      return f, -1
  end
  return base_case
end

return _recursive_reduction(f, g, A=A, reps=reps, base_case=base_case, open=open)
end

#=
n = 5
G1 = g1_named(n)
u = magic_unitary(n)
g = G1.gs
bg934 = bg(9,3,4,u=u)
normal_form_with_rep_choices(bg934, g)
recursive_reduction(bg934,g)


=#
function normal_form_with_rep_choices(
  f::Generic.FreeAssociativeAlgebraElem{T},
  g::Vector{Generic.FreeAssociativeAlgebraElem{T}},
  len::Int=-1
) where {T}

  A = parent(f)
  reps = Tuple{elem_type(base_ring(A)),Generic.FreeAssociativeAlgebraElem{T},Generic.FreeAssociativeAlgebraElem{T},Generic.FreeAssociativeAlgebraElem{T}}[]
  s = length(g)
  @label first
  i = 1
  @label again
  iszero(f) && return f, reps
  ok, ml, mr = Generic.word_divides_leftmost(f.exps[1], g[i].exps[1])
  #Figure out all choices
  choices = get_all_divisors(f, g)
  if length(choices) > 1
    @info "After $(length(reps)) steps, there are $(length(choices)) choices."
  end

  if !ok && i < s
    i += 1
    @goto again
  end

  if ok && !iszero(f)
    qi = divexact(f.coeffs[1], g[i].coeffs[1])
    new_f = Generic._sub_rest(f, Generic.mul_term(qi, ml, g[i], mr), 1) # enforce lt cancelation

    mul_l = one(A)
    mul_r = one(A)

    if length(ml) > 0
      mul_l = prod([A[i] for i in ml])
    end
    if length(mr) > 0
      mul_r = prod([A[i] for i in mr])
    end
    push!(reps, (qi, mul_l, g[i], mul_r))

    #break if computation limit is set
    if len>=0 && length(reps)>=len
      @warn "Not a full reduction! Only $len steps were computed."
      return f,reps
    end

    f = new_f
    #Check if the leading monomial has been killed

    @goto first
  else
    return f, reps
  end
end


