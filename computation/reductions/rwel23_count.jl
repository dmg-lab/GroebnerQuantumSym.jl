#= The counting Argument
The Normal Form

Syntax ∑_{i1,i2,i3,i4} u(i1,i2,i3,i4)

s1 := 2:n × {2} × 2:n × 3:n
filter: i3 ≠ i1
coeff := 1

d1 := 2:n × 3:n × 2:n × {1}
filter: i3 ≠ i1
coeff := -1

d2 := {2} × 2:n × 3:n × (4:n ∪ {2})
filter: i4 ≠ i2
coeff := -1

s2 := 3:n × 2:n × {1} × (4:n ∪ {2})
filter: i4 ≠ i2
coeff := 1

d3 := 3:n × {2} × 1:n × 2:n
filter: i4 ≠ i2, i3 ≠ i1
coeff := -1

s4 := {2} × 3:n × 2:n × 1:n
filter: i4 ≠ i2, i3 ≠ i1
coeff := 1

s5 := 3:n × 3:n × 2:n × 1:n
filter: i4 ≠ i2, i3 ≠ i1
coeff:= 1

d6 := 3:n × 3:n × 1:n × 2:n
filter: i4 ≠ i2, i3 ≠ i1
coeff:= -1

d7 := {2} × 4:n × 3:n × {3}
filter: none
coeff:= -1

s7 := 3:n × 4:n × {1} × {3}
filter: none
coeff:= 1

d0 := {2} × {2} × 3:n × {3}
filter: none
coeff:= -1

s0 := 3:n × {2} × {1} × {3}
filter: none
coeff:= 1

=#

ElementType = Union{Symbol, Int}
TaggedTuple = Tuple{NTuple{N,ElementType}, Vector{Symbol}} where N

s1 = Iterators.product([2,3,:g], [2], [2,3,:g], [3,:g])
s1 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], s1)
s1tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, Symbol[:g31]) : (x,Symbol[]) for x in s1]

d1 = Iterators.product([2,3,:g], [3,:g], [2,3,:g], [1])
d1 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], d1)
d1tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, Symbol[:g31]) : (x,Symbol[]) for x in d1]

d2 = Iterators.product([2], [2,3,:g], [3,:g], [2,:g])
d2 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], d2)
d2tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in d2]

s2 = Iterators.product([3,:g], [2,3,:g], [1], [2,:g])
s2 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], s2)
s2tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in s2]

d3 = Iterators.product([3,:g], [2], [1,2,3,:g], [2,3,:g])
d3 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], d3)
d3 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], d3)
d3tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in d3]
d3tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, push!(y, :g31)) : (x, y) for (x,y) in d3tagged]

s4 = Iterators.product([2], [3,:g], [2,3,:g], [1,2,3,:g])
s4 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], s4)
s4 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], s4)
s4tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in s4]
s4tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, push!(y, :g31)) : (x, y) for (x,y) in s4tagged]

s5 = Iterators.product([3,:g], [3,:g], [2,3,:g], [1,2,3,:g])
s5 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], s5)
s5 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], s5)
s5tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in s5]
s5tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, push!(y, :g31)) : (x, y) for (x,y) in s5tagged]

d6 = Iterators.product([3,:g], [3,:g], [1,2,3,:g], [2,3,:g])
d6 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], d6)
d6 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], d6)
d6tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in d6]
d6tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, push!(y, :g31)) : (x, y) for (x,y) in d6tagged]


d7 = Iterators.product([2], [:g], [3,:g], [3])
d7tagged = TaggedTuple[(x, Symbol[]) for x in d7]

s7 = Iterators.product([3,:g], [:g], [1], [3])
s7tagged = TaggedTuple[(x, Symbol[]) for x in s7]

d0 = Iterators.product([2], [2], [3,:g], [3])
d0tagged = TaggedTuple[(x, Symbol[]) for x in d0]

s0 = Iterators.product([3,:g], [2], [1], [3])
s0tagged = TaggedTuple[(x, Symbol[]) for x in s0]

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

all = [
    (s1tagged, 1),
    (d1tagged, -1),
    (d2tagged, -1),
    (s2tagged, 1),
    (d3tagged, -1),
    (s4tagged, 1),
    (s5tagged, 1),
    (d6tagged, -1),
    (d7tagged, -1),
    (s7tagged, 1),
    (d0tagged, -1),
    (s0tagged, 1)
]

count_cases(all;debug=false)

