
using Oscar.AbstractAlgebra.Generic




"""
    normal_form_noRedTail(f::Generic.FreeAssAlgElem{T}, g::Vector{Generic.FreeAssAlgElem{T}}) where T<:FieldElement

Compute the normal form of `f` with respect to the leading term of the elements of `g`. Without tail reduction.

```jldoctest
using Oscar
R, (x, y) = free_associative_algebra(QQ, ["x", "y"]);
f = 2x^2 + 3x*y + 4y^2;
g = [x^2 + y^2];

normal_form_noRedTail(f, g)

# output

3*x*y + 2*y^2
```

"""
function normal_form_noRedTail(
    f::Generic.FreeAssAlgElem{T},
    g::Vector{Generic.FreeAssAlgElem{T}},
) where T<:FieldElement
    otp = "="^8*"Normal Form without Tail Reduction"*"="^8 * "\n"
    otp *= "f = "*string(f)*"\n"
    otp *="-"^50 * "\n"
    s = length(g)
    @label first
    i = 1
    @label again 
    iszero(f) && return f 
    ok, ml, mr = Generic.word_divides_leftmost(f.exps[1], g[i].exps[1])
    if !ok && i < s
        i += 1
        @goto again
    end

    if ok && !iszero(f)
        qi = divexact(f.coeffs[1], g[i].coeffs[1])
        f = Generic._sub_rest(f, Generic.mul_term(qi, ml, g[i], mr), 1) # enforce lt cancelation
        otp *= "Reducing f by "*string(qi)*"*"*string(g[i])*"\n"
        otp *= "new f = "*string(f)*"\n"
        otp *="-"^50 * "\n"

        #Check if the leading monomial has been killed
        
        @goto first
    else
        otp *= "="^4*"End of Normal Form without Tail Reduction"*"="^5 * "\n"
        println(otp)
        return f
    end
end

"""
    normal_form_with_rep(f::Generic.FreeAssAlgElem{T}, g::Vector{Generic.FreeAssAlgElem{T}}) where T<:FieldElement

Compute the normal form of `f` with respect to the leading term of the elements of `g`. Without tail reduction. 
It stores the reduction steps in a dictionary.

```jldoctest
using Oscar
R, (x, y) = free_associative_algebra(QQ, ["x", "y"]);
f = 2x^2 + 3x*y*x^2 + 4y^2;
g = [x^2, y^2, x*y];

zr,  dct = normal_form_with_rep(f, g);
dct
zr
# output

Dict{AbstractAlgebra.Generic.FreeAssAlgElem, AbstractAlgebra.Generic.FreeAssAlgElem{QQFi
eldElem}} with 2 entries:
  x^2 => 2*x^2
  y^2 => 4*y^2
```

"""
function normal_form_with_rep(
    f::Generic.FreeAssAlgElem{T},
    g::Vector{Generic.FreeAssAlgElem{T}},
) where {T}

    rep_dict = Vector{Tuple{Generic.FreeAssAlgElem, Generic.FreeAssAlgElem{T}}}()
    s = length(g)
    @label first
    i = 1
    @label again 
    iszero(f) && return f , rep_dict
    ok, ml, mr = Generic.word_divides_leftmost(f.exps[1], g[i].exps[1])
    if !ok && i < s
        i += 1
        @goto again
    end

    if ok && !iszero(f)
        qi = divexact(f.coeffs[1], g[i].coeffs[1])
        new_f = Generic._sub_rest(f, Generic.mul_term(qi, ml, g[i], mr), 1) # enforce lt cancelation

        push!(rep_dict, (g[i], f-new_f))
        println("new_f: $new_f")
        
        f = new_f
        #Check if the leading monomial has been killed
        
        @goto first
    else
        return f , rep_dict
    end
end

"""
    interreduce!_noRedTail(g::Vector{Generic.FreeAssAlgElem{T}}) where T

Interreduce a list of elements without tail reduction.

```
using Oscar
using QuantumAutomorphismGroups
relat , _ = getQuantumPermutationGroup(4);

red_relat = interreduce!_noRedTail(relat)
length(red_relat)

# output

63
```

"""
function interreduce!_noRedTail(g::Vector{Generic.FreeAssAlgElem{T}}) where T<:FieldElement
    i = 1
    while length(g) > 1 && length(g) >= i
        r = normal_form_noRedTail(g[i], g[1:end .!= i])
        if iszero(r)
            deleteat!(g, i)
        elseif g[i] != r
            g[i] = r
            i = 1
        else
            i += 1
        end
    end
    return g
end

"""
    get_obstruction_pairs(g_old = Vector{Generic.FreeAssAlgElem{T}}) where T<:FieldElement

Compute the obstruction pairs of a list of elements.

```jldoctest
using Oscar
using QuantumAutomorphismGroups
relat , _ = getQuantumPermutationGroup(4);
get_obstruction_pairs(relat);
first_poly, second_poly, Spoly_red, _ = get_obstruction_pairs(relat)[1];
Spoly_red

# output

u[4,2]*u[3,3] + u[4,2]*u[3,4] + u[4,3]*u[3,2] + u[4,3]*u[3,3] + u[4,3]*u[3,4] + u[4,4]*u[3,2] + u[4,4]*u[3,3] + u[4,4]*u[3,4] - u[3,2] - u[3,3] - u[3,4] + u[4,1]
```

"""
function get_obstruction_pairs(g_old::Vector{Generic.FreeAssAlgElem{T}}) where T<:FieldElement
    g = copy(g_old)
    # step 1 from Thm. 5.2.12 Noncommutative Groebner Bases and Applications, Xingqiang Xiu
    obstruction_queue = Generic.get_obstructions(g)
    ans=[]
    while !isempty(obstruction_queue) 
        obstruction = popfirst!(obstruction_queue)[1]
        # step3
        S = Generic.s_polynomial(obstruction)
        Sp = normal_form_noRedTail(S, g) # or normal_form_weak
        if iszero(Sp)
            continue
        else
            push!(ans,(obstruction.first_poly,obstruction.second_poly,Sp,S)) 
        end
    end
    return ans
end


"""
    getQuantumRelationsByType(n::Int)

Compute the quantum relations of a quantum permutation group of degree `n` by type.

```julia
using Oscar
cs, rs, id, wel, inj = getQuantumRelationsByType(4)

```
"""
function getQuantumRelationsByType(n::Int)
    relat , relat_by_type, u, _ = getQuantumPermutationGroup(n);
    cs = typeof(relat[1])[]
    for ele in relat_by_type[:col_sum]
        push!(cs, ele)
    end
    rs = typeof(relat[1])[]
    for ele in relat_by_type[:row_sum]
        push!(rs, ele)
    end
    id = permutedims(reshape(relat_by_type[:idempotent],(size(u)[1],size(u)[2])),[2,1])

    #Define wels to be vector of even index entries of relat_by_type[:zero_divisor] inj
    wels = typeof(relat[1])[]
    inj = typeof(relat[1])[]
    for i in 1:1:length(relat_by_type[:zero_divisor])
        if i % 2 == 0
            push!(wels,relat_by_type[:zero_divisor][i])
        else
            push!(inj,relat_by_type[:zero_divisor][i])
        end
    end
    reshaped_wels = Array{typeof(relat[1]),3}(undef, (size(u)[1],size(u)[1],size(u)[2]))

    sk = div(length(wels),n)
    for h in 0:length(wels)-1
        j = div(h,sk)+1
        h1 = rem(h,sk)+1
        i = div(h1-1,n-1) +1
        k = rem(h1-1,n-1) +1
        if k >= i
            k += 1
        end
        reshaped_wels[i,j,k] = wels[h+1]
    end

    reshaped_inj = Array{typeof(relat[1]),3}(undef, (size(u)[1],size(u)[1],size(u)[2]))

    sk = div(length(inj),n)
    for h in 0:length(inj)-1
        j = div(h,sk)+1
        h1 = rem(h,sk)+1
        i = div(h1-1,n-1) +1
        k = rem(h1-1,n-1) +1
        if k >= i
            k += 1
        end
        reshaped_inj[j,i,k] = inj[h+1]
    end
    return cs, rs, id, reshaped_inj, reshaped_wels
end




