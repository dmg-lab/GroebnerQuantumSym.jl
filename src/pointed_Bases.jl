#Pointed Bases
using Oscar
using QuantumAutomorphismGroups
using Oscar.AbstractAlgebra.Generic

export 
    PointedSet,
    getPointedSet,
    parent,
    representative,
    getPointedSetRelation,
    getPointedSetRelations




mutable struct PointedSet
    point::Int
    base::Vector{Int}
    rep::Generic.FreeAssociativeAlgebraElem{T} where T<:FieldElem
    function PointedSet(point::Int,base::Vector{Int}, rep::Generic.FreeAssociativeAlgebraElem{T}) where T<:FieldElem
        @assert point in base
        new(point,base,rep)
    end
end

QuantumGB.parent(p::PointedSet) = Oscar.parent(p.rep)
QuantumGB.representative(p::PointedSet) = p.rep

#=
using Oscar
using QuantumAutomorphismGroups

M = uniform_matroid(2,4)
ps = getPointedSet(M)


getPointedSetRelations(M)

=#
function getPointedSet(M::Matroid;structure::Symbol=:bases)
    sets  = eval(structure)(M)
    r = rank(M)
    relat, relat_, u, A = getQuantumPermutationGroup(r * length(sets))

    answer = PointedSet[]
    j = 1
    for base in sets
        pb = PointedSet[]
        for i in base
            push!(pb,PointedSet(i,base,A[j]))
            j += 1
        end
        append!(answer,pb)
    end
    return answer, (relat, relat_, u, A)
end

function getPointedSetRelation(p1::PointedSet,p2::PointedSet)
    if p1.point == p2.point
        if p1.base == p2.base
            return 2
        else
            return 0
        end
    else
        if p1.point in p2.base && p2.point in p1.base
            return 1
        else
            return -1
        end 
    end
end

function getPointedSetRelations(M::Matroid,structure::Symbol=:bases)
    sets, (relat,_,_,_)  = getPointedSet(M,structure=structure)


    return relat
end

B1 = [1,2,3,5]
B2 = [2,3,4,5]

# Get the Vector of all 4 elment subsets of [1,2,3,4,5]

function subsets(set::Vector{Vector{Int}})
    n = length(set)
    result = typeof(set)[]
    for i in 0:(2^n - 1)
        subset = []
        for j in 1:n
            if ((i >> (j - 1)) & 1) != 0
                push!(subset, set[j])
            end
        end
        push!(result, subset)
    end
    return result
end

function subsets_of_size(set::Vector{Int},size::Int)
    n = length(set)
    result = typeof(set)[]
    for i in 0:(2^n - 1)
        subset = []
        for j in 1:n
            if ((i >> (j - 1)) & 1) != 0
                push!(subset, set[j])
            end
        end
        if length(subset) == size
            push!(result, subset)
        end
    end
    return result
end

subsets_of_size(set::Int,size::Int) = subsets_of_size(collect(1:set),size)

sets = subsets_of_size(7,4)
all_subsets = subsets(sets)

for s in all_subsets
    if length(s) ==  0 
        continue
    end
    matroid_from_bases(s,5,check=true)
end


M = matroid_from_bases([B1,B2],5)










