
export run_tests, count_cases, TaggedTuple, ElementType


ElementType = Union{Symbol, Int}
TaggedTuple = Tuple{NTuple{N,ElementType}, Vector{Symbol}} where N

function count_cases(cases::Array{TaggedTuple}, dct::Dict{TaggedTuple,Int}=Dict(),sgn::Int=1;debug=false)
  for c in cases
    if haskey(dct, c)
      dct[c] += sgn
      if dct[c] == 0 && !debug
        delete!(dct, c)
      end
    else
      dct[c] = sgn
    end
  end
  return dct
end

function count_cases(cases::Vector{Tuple{Array{TaggedTuple},Int}}, dct::Dict{TaggedTuple,Int}=Dict{TaggedTuple,Int}();debug=false)
  for (c, sgn) in cases
    dct = count_cases(c, dct, sgn;debug=debug)
  end
  return dct
end

function count_cases(cases::Vector{Tuple{Vector{TaggedTuple},Int}}, dct::Dict{TaggedTuple,Int}=Dict{TaggedTuple,Int}();debug=false)
  for (c, sgn) in cases
    dct = count_cases(c, dct, sgn;debug=debug)
  end
  return dct
end

function count_cases(cases::Vector{Tuple{Matrix{TaggedTuple},Int}}, dct::Dict{TaggedTuple,Int}=Dict{TaggedTuple,Int}();debug=false)
  for (c, sgn) in cases
    for x in c
      dct = count_cases(x, dct, sgn;debug=debug)
    end
  end
  return dct
end

function run_tests()
  include("../test/runtests.jl");
  return nothing
end
