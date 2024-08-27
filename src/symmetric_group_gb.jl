#Test

export rwel, rinj, magic_unitary, inj, wel, normed, g0



function magic_unitary(n::Int=4)
  A, u = free_associative_algebra(QQ, :u => (1:n, 1:n))
  return u
end

Base.length(A::FreeAssAlgebra) = length(A.S)

magic_unitary(A::FreeAssAlgebra) = magic_unitary(Int(sqrt(length(A))))
magic_unitary(f::Generic.FreeAssAlgElem) = magic_unitary(parent(f))


function rwel(k::Int, j::Int, h::Int=2, v::Int=2; u::Matrix{Generic.FreeAssAlgElem{T}}=magic_unitary()) where T
  n = size(u)[1]
  return sum([u[2, k] * u[s, j] for s in h+1:n]) - sum([u[s, k] * u[1, j] for s in v+1:n]) + u[1, j] - u[2, k]
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


function g0(n::Int=4, check::Bool=false)
  rel, relat_by_type, u, A = getQuantumPermutationGroup(n)

  cs = relat_by_type[:col_sum]
  rs = relat_by_type[:row_sum][2:end]

  #Idempotent relations but not the ones that contain i,j = 1
  ip = relat_by_type[:idempotent]
  nonip = filter(x -> lm(x) in vcat([u[1, j]^2 for j in 1:n], [u[i, 1]^2 for i in 1:n]), ip)
  ip_f = filter(x -> !(x in nonip), ip)

  wels = [u[i, j] * u[i, k] for i in 2:n, j in 2:n, k in 2:n if j != k]
  inj = [u[i, j] * u[k, j] for i in 2:n, j in 2:n, k in 2:n if i != k]

  rwels = reshape([rwel(k, j; u = u) for j in 2:n, k in 2:n if j != k], (n - 1) * (n - 2))
  rinjs = reshape([rinj(k, j; u = u) for j in 2:n, k in 2:n if j != k], (n - 1) * (n - 2))

  #filter out shit
  rwels_min = filter(x -> lm(x) != u[2, 2] * u[3, 3], rwels)

  gens = vcat(ip_f, rs, cs, wels, inj, rwels_min, rinjs)

  return gens
end






#=
using Oscar
G0 = g0(4)

gb1 = groebner_basis(G0)
normed.(gb!)
rel, _ , _ , _ = getQuantumPermutationGroup(4)
gb2 = groebner_basis(rel)


groebner_basis(gb1)
AbstractAlgebra.interreduce!(gb1)


max_degree(gb!)


=#




