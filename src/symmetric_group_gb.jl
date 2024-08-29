#Test

export rwel,
  rinj,
  magic_unitary,
  inj,
  wel,
  ip, 
  normed,
  g0,
  indexed_name,
  row_sum,
  col_sum,
  gg,
  g0_extended,
  red_row_sum,
  red_col_sum,
  red_sums,
  red_guy,
  rinj_col_sum,
  rwel_col_sum, 
  bg1,
  bg2,
  bg3,
  bg4,
  bg5,
  bg6,
  bg7,
  bg8,
  bg9,
  bg10,
  bg11,
  bg12,
  bg13,
  bg14



function magic_unitary(n::Int=4)
  A, u = free_associative_algebra(QQ, :u => (1:n, 1:n))
  return u
end

Base.length(A::FreeAssAlgebra) = length(A.S)

function magic_unitary(A::FreeAssAlgebra)
  n = Int(sqrt(length(A)))
  return permutedims(reshape([A[i] for i in 1:length(A)], (n, n)), [2, 1])
end

magic_unitary(f::Generic.FreeAssAlgElem) = magic_unitary(parent(f))
function magic_unitary(v::Vector{Generic.FreeAssAlgElem{T}}) where T
  @assert all([parent(x) == parent(v[1]) for x in v]) "All elements must be in the same algebra"
  @assert length(v) > 0 "Empty vector"
  return magic_unitary(parent(v[1])) 
end



function rwel(k::Int, j::Int, h::Int=2, v::Int=2; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  n = size(u)[1]
  return sum([u[2, k] * u[s, j] for s in h+1:n];init=zero(T)) - sum([u[s, k] * u[1, j] for s in v+1:n]) + u[1, j] - u[2, k]
end

function rinj(k::Int, j::Int,  h::Int=2, v::Int=2; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  n = size(u)[1]
  return sum([u[k, 2] * u[j, s] for s in h+1:n]) - sum([u[k, s] * u[j, 1] for s in v+1:n]) + u[j, 1] - u[k, 2]
end

function wel(i::Int, j::Int, k::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  @assert 1 <= i <= size(u)[1] && 1 <= j <= size(u)[2] && 1 <= k <= size(u)[2] "Indices out of bounds"
  return u[i, j] * u[i, k]
end

function inj(i::Int, j::Int, k::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  @assert 1 <= i <= size(u)[1] && 1 <= j <= size(u)[2] && 1 <= k <= size(u)[1] "Indices out of bounds"
  return u[i, j] * u[k, j]
end

function ip(i::Int, j::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  @assert 1 <= i <= size(u)[1] && 1 <= j <= size(u)[2] "Indices out of bounds"
  return u[i, j] * u[i, j] - u[i, j]
end


#=
u = magic_unitary()
f = -4 * u[1,2] * u[2,3] + u[1,2]
normed(f)
=#
function normed(f::Generic.FreeAssAlgElem)
  c1 = f.coeffs[1]
  T = typeof(c1)
  normed_c1 = T(sqrt(c1^2))
  f.coeffs[1] = normed_c1
  return f
end


function g0(n::Int=4; names=false)
  rel, relat_by_type, u, A = getQuantumPermutationGroup(n)

  cs = relat_by_type[:col_sum]
  cs_names = indexed_name("cs", 1:length(cs))
  rs = relat_by_type[:row_sum][2:end]
  rs_names = indexed_name("rs", 2:length(rs)+1)

  #Idempotent relations but not the ones that contain i,j = 1
  ip = reshape([u[i, j]^2 - u[i, j] for i in 2:n, j in 2:n], (n-1)^2)
  ip_names = indexed_name("ip", reshape([parse(Int,"$i$j") for i in 2:n, j in 2:n], (n-1)^2))
  wels = [u[i, j] * u[i, k] for i in 2:n, j in 2:n, k in 2:n if j != k]
  wels_names = indexed_name("wel", [parse(Int,"$i$j$k") for i in 2:n, j in 2:n, k in 2:n if j != k])
  inj = [u[i, j] * u[k, j] for i in 2:n, j in 2:n, k in 2:n if i != k]
  inj_names = indexed_name("inj", [parse(Int,"$i$j$k") for i in 2:n, j in 2:n, k in 2:n if i != k])
  rwels = reshape([rwel(k, j; u = u) for j in 2:n, k in 2:n if j != k], (n - 1) * (n - 2))
  rwels_names = indexed_name("rwel", reshape([parse(Int,"$j$k") for j in 2:n, k in 2:n if j != k], (n - 1) * (n - 2)))
  rinjs = reshape([rinj(k, j; u = u) for j in 2:n, k in 2:n if j != k], (n - 1) * (n - 2))
  rinjs_names = indexed_name("rinj", reshape([parse(Int,"$j$k") for j in 2:n, k in 2:n if j != k], (n - 1) * (n - 2)))
  filter!(x -> x != "rwel₂₃", rwels_names)
  #filter out shit
  rwels_min = filter(x -> lm(x) != u[2, 2] * u[3, 3], rwels)

  gens = vcat(ip, rs, cs, wels, inj, rwels_min, rinjs)
  gens_names = vcat(ip_names,rs_names,cs_names,wels_names,inj_names, rwels_names,rinjs_names)

  if names
    @assert length(gens) > 0
    @assert length(gens) == length(gens_names)
    names_dct = Dict{typeof(gens[1]),String}()
    for i in eachindex(gens)
      names_dct[gens[i]] = gens_names[i]
    end
    return gens, names_dct
  end

  return gens
end

row_sum(i::Int, u::Matrix{Generic.FreeAssAlgElem{T}} where T =magic_unitary()) = sum([u[i,x] for x in 1:size(u)[1]]) - one(parent(u[1,1]))
col_sum(i::Int, u::Matrix{Generic.FreeAssAlgElem{T}} where T =magic_unitary()) = sum([u[x,i] for x in 1:size(u)[1]]) - one(parent(u[1,1]))


function gg(x::Int, y::Int, u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  n = size(u)[1]
  gg_base =  u[x,y] * (sum([row_sum(i,u) for i in 2:n if i != x]) - sum([col_sum(i,u) for i in 2:n if i != y]))
  gg_plus = sum([wel(x,y,i; u = u) for i in 2:n if i != y]) - sum([inj(x,y,i; u = u) for i in 2:n if i != x])
  return gg_base + gg_plus
end

function indexed_name(name::String, numbers::Matrix{Int})
  numbers = vcat(numbers...)
  return indexed_name(name,numbers) 
end

function indexed_name(name::String, numbers::UnitRange{Int})
  indexed_name(name,collect(numbers))
end

function indexed_name(name::String, numbers::Vector{Int})
  return [name * join(["₀₁₂₃₄₅₆₇₈₉"[3*i+1] for i in reverse(digits(n))], "") for n in numbers]
end

function red_row_sum(x::Int, y::Int, k::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return u[x,y] * row_sum(k,u) - inj(x,y,k;u)
end
function red_col_sum(x::Int, y::Int, k::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return u[x,y] * col_sum(k,u) - wel(x,y,k;u)
end

function red_sums(u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  n = size(u)[1]
  rrs = [sum([red_row_sum(i,j,k;u) for k in 2:n if i != k]) for i in 2:n for j in 2:n]
  rcs = [sum([red_col_sum(i,j,k;u) for  k in 2:n if j != k]) for i in 2:n for j in 2:n]
  return vcat(rrs,rcs)
end

function red_guy(x::Int, y::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  n = size(u)[1]
  return sum([rinj(k,y;u) for k in 2:n if k != y]) - sum([red_row_sum(y,x,k;u) for k in 2:n if k != y])
end

function rinj_col_sum(x::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  n = size(u)[1]
  return sum([rinj(k,x;u) for k in 2:n if k != x])
end

function rwel_col_sum(x::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  n = size(u)[1]
  return sum([rwel(k,x;u) for k in 2:n if k != x])
end



function g0_extended(n::Int; names=false)
  if names
    gns, gns_names = g0(n, names=true)
    u = magic_unitary(gns)
    ggs = [gg(i,j,u) for i in 2:n for j in 2:n]
    ggs_names = indexed_name("gg", reshape([parse(Int,"$i$j") for i in 2:n, j in 2:n], (n-1)^2))

    red_sums_v = red_sums(u)
    red_sums_names = indexed_name("rrs", [parse(Int,"$i$j") for i in 2:n for j in 2:n])
    red_cols_names = indexed_name("rcs", [parse(Int,"$i$j") for i in 2:n for j in 2:n])
    red_sums_names = vcat(red_sums_names, red_cols_names)

    rinj_wel_col_sums = vcat([rinj_col_sum(i;u) for i in 2:n], [rwel_col_sum(i;u) for i in 4:n])
    rinj_wel_col_sums_names = vcat(indexed_name("rinjcs", 2:n), indexed_name("rwelcs", 4:n))
    for i in eachindex(rinj_wel_col_sums_names)
      gns_names[rinj_wel_col_sums[i]] = rinj_wel_col_sums_names[i]
    end
    
    #rgs = [red_guy(i,j;u) for i in 2:n for j in 2:n if i != j]
    #rgs_names = indexed_name("mrrs", [parse(Int,"$i$j") for i in 2:n for j in 2:n if i != j])
    #for i in eachindex(rgs_names)
    #  gns_names[rgs[i]] = rgs_names[i]
    #end


    gns = vcat(ggs,rinj_wel_col_sums,red_sums_v, gns)


    for i in eachindex(ggs_names)
      gns_names[ggs[i]] = ggs_names[i]
    end
    for i in eachindex(red_sums_names)
      gns_names[red_sums_v[i]] = red_sums_names[i]
    end
    return gns, gns_names
  else
    gns = g0(n)
    u = magic_unitary(gns)
    ggs = reshape([gg(i,j,u) for i in 2:n, j in 2:n], n^2)
    red_sums_v = red_sums(u)
    gns = vcat(ggs,red_sums_v, gns)
    return gns
  end
end

#=
using Oscar
G0, names = g0(4,names=true)
u = magic_unitary(G0)
bg1 = rwel(2,3; u = u)
normal_form(bg1,G0)
names

G0, names = g0_extended(6,names=true)
u = magic_unitary(G0)
bg1 = rwel(2,3; u = u)
normal_form(bg1,G0)
names

mrrs = mod_rrs(2,3; u = u)

r, v = normal_form_with_rep(mrrs,G0);
r != 0 && error("r != 0");
x1 = reps_vector_to_poly(v,names);
print(x1)

digits(5)
indexed_name("test", [2,4,5,15])

groebner_basis(gb1)
AbstractAlgebra.interreduce!(gb1)


max_degree(gb!)


=#

function bg1(i::Int, j::Int, k::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return inj(i,2,k,u=u)*u[j, 3] -  u[i, 2]*rinj(k,j,u=u)
end

function bg2(i::Int, j::Int, k::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return u[k,2]*inj(j,3,i,u=u) -  rinj(k,j,u=u)*u[i,3]
end

function bg3(k::Int, j::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return ip(2,k,u=u)*u[3,j] -  u[2,k]*rwel(k,j,u=u)
end

function bg4(k::Int, j::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return u[2,k]*ip(3,j,u=u) -  rwel(k,j,u=u)*u[3,j]
end

function bg5(k::Int, j::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return ip(k,2,u=u)*u[j,3] -  u[k,2]*rinj(k,j,u=u)
end

function bg6(k::Int, j::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return u[k,2]*ip(k,2,u=u) -  rinj(k,j,u=u)*u[j,3]
end

function bg7(j::Int, k::Int, h::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return wel(2,j,k,u=u)*u[3, h] -  u[2, j]*rwel(k,h,u=u)
end

function bg8(k::Int, j::Int, h::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return u[2,k]*wel(3,j,h,u=u) -  rwel(k,j,u=u)*u[3,h]
end

function bg9(k::Int, j::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return rinj(k,2,u=u)*u[3,j] -  u[k,2]*rwel(3,j,u=u)
end

function bg10(j::Int, k::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return u[2,j]*rinj(3,k,u=u) -  rwel(j,2,u=u)*u[k,3]
end

function bg11(k::Int, i::Int, j::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return inj(k,i,2,u=u)*u[3,j] -  u[k,i]*rwel(i,j,u=u)
end

function bg12(k::Int, j::Int, h::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return u[2,k]*inj(3,i,j,u=u) -  rwel(k,i,u=u)*u[j,i]
end

function bg13(i::Int, j::Int, k::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return wel(i,j,2,u=u)*u[k, 3] -  u[i,j]*rinj(i,k,u=u)
end

function bg14(k::Int, i::Int, j::Int; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  return u[k,2]*wel(i,3,j,u=u) -  rinj(k,i,u=u)*u[i,j]
end
