#Checking the conjecture
export lm_gextra, lm_gb, lm_conjecture, check_conjecture


#=


=#
function lm_gextra(u::Matrix{Generic.FreeAssociativeAlgebraElem{T}}=magic_unitary()) where T<:FieldElem
  n = size(u)[1]
  e1 = [u[2,i]*u[4,j]*u[3,k] for i=3:n for j=3:n for k=2:n if i!=j && j!=k]
  e1_s = [u[i,2]*u[j,4]*u[k,3] for i=3:n for j=3:n for k=2:n if i!=j && j!=k]
  e14 = [u[2,i]*u[4,2]*u[3,k] for i=3:n for k=4:n]
  e14_s = [u[i,2]*u[2,4]*u[k,3] for i=3:n for k=4:n]

  e2 = [u[2,2]*u[4,j]*u[3,k] for j in 5:n for k in 2:n if j!=k]
  e2_s = [u[2,2]*u[j,4]*u[k,3] for j in 5:n for k in 2:n if j!=k]

  e3 = [u[2,2]*u[4,4]*u[3,k] for k in 2:n if 4!=k && 3!=k]
  e3_s = [u[2,2]*u[4,4]*u[k,3] for k in 2:n if 4!=k && 3!=k]

  e4 = [u[2,2]*u[4,4]*u[3,3]]
  E = vcat(e1,e1_s,e14,e14_s,e2,e2_s,e3,e3_s,e4)
  @assert length(E) == gb_count(n) - g0_count(n) "Something went wrong: $(length(E)) != $(gb_count(n) - g0_count(n))"
  return sort!(E)
end

lm_gextra(n::Int) = lm_gextra(magic_unitary(n))

function lm_gb(n::Int=4;deg_bound::Int=-1, return_magic_unitary::Bool=false)
  G0 = g0(n)
  gb = Oscar.groebner_basis(G0, deg_bound)

  lms_gb = sort!(Oscar.leading_monomial.(gb))
  if return_magic_unitary
    return lms_gb, magic_unitary(gb)
  end
  return lms_gb
end


lm_conjecture(n::Int=4) = lm_conjecture(magic_unitary(n))

function lm_conjecture(u::Matrix{Generic.FreeAssociativeAlgebraElem{T}}) where T<:FieldElem
  n = size(u)[1]
  G0 = g0(n)
  lme = lm_gextra(u)
  lmg0 = Oscar.leading_monomial.(G0)
  ret = sort!(vcat(lme,lmg0)) 
  @assert length(ret) == gb_count(n) "Something went wrong: $(length(ret)) != $(gb_count(n))"
  return ret
end

function check_conjecture(n::Int=4; deg_bound::Int=-1)
  @assert n > 3 "n must be greater than 3"
  G0 = g0(n)
  gb = Oscar.groebner_basis(G0, deg_bound)
  u = magic_unitary(gb)
  lms_gb = sort!(Oscar.leading_monomial.(gb))
  conj = lm_conjecture(u)

  otp = "" 

  length(gb) == length(conj) || (otp *= "Lengths of GB and Conjecture are different\n") 
  isempty(setdiff(lms_gb,conj)) || (otp *= "GB is not contained in the conjecture\n")
  isempty(setdiff(conj,lms_gb)) || (otp *= "Conjecture is not contained in the GB\n")
  if isempty(otp) 
      otp *=  "The conjecture is correct for n = $n: length(gb) = $(length(lms_gb)) = length(conj)\n"
  end
  otp *= "The following groebner basis was computed(n=$n):\n"
  for ele in gb
    otp *= "$ele\n"
  end

  return otp
end


#=
x,u = lm_gb(4,return_magic_unitary=true)
y = lm_conjecture(u)
setdiff(x,y)

print(check_conjecture(4))

=#





