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
@testset "rwel23_old" begin
println("""
Testing the general case of rwel23 (old version)...""")


@testset "Degree 2" begin

s1 = Iterators.product([2,3,:g], [2], [2,3,:g], [3,:g]);
s1 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], s1);
s1tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, Symbol[:g31]) : (x,Symbol[]) for x in s1];

d1 = Iterators.product([2,3,:g], [3,:g], [2,3,:g], [1]);
d1 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], d1);
d1tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, Symbol[:g31]) : (x,Symbol[]) for x in d1];

d2 = Iterators.product([2], [2,3,:g], [3,:g], [2,:g]);
d2 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], d2);
d2tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in d2];

s2 = Iterators.product([3,:g], [2,3,:g], [1], [2,:g]);
s2 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], s2);
s2tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in s2];

d3 = Iterators.product([3,:g], [2], [1,2,3,:g], [2,3,:g]);
d3 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], d3);
d3 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], d3);
d3tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in d3];
d3tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, push!(y, :g31)) : (x, y) for (x,y) in d3tagged];

s4 = Iterators.product([2], [3,:g], [2,3,:g], [1,2,3,:g]);
s4 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], s4);
s4 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], s4);
s4tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in s4];
s4tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, push!(y, :g31)) : (x, y) for (x,y) in s4tagged];

s5 = Iterators.product([3,:g], [3,:g], [2,3,:g], [1,2,3,:g]);
s5 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], s5);
s5 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], s5);
s5tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in s5];
s5tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, push!(y, :g31)) : (x, y) for (x,y) in s5tagged];

d6 = Iterators.product([3,:g], [3,:g], [1,2,3,:g], [2,3,:g]);
d6 = Iterators.filter(x -> x[4] == :g || x[2] == :g || x[4] != x[2], d6);
d6 = Iterators.filter(x -> x[3] == :g || x[1] == :g || x[3] != x[1], d6);
d6tagged = TaggedTuple[x[4] == :g && x[2] == :g ? (x, Symbol[:g42]) : (x,Symbol[]) for x in d6];
d6tagged = TaggedTuple[x[3] == :g && x[1] == :g ? (x, push!(y, :g31)) : (x, y) for (x,y) in d6tagged];

d7 = Iterators.product([2], [:g], [3,:g], [3]);
d7tagged = TaggedTuple[(x, Symbol[]) for x in d7];

s7 = Iterators.product([3,:g], [:g], [1], [3]);
s7tagged = TaggedTuple[(x, Symbol[]) for x in s7];

d0 = Iterators.product([2], [2], [3,:g], [3]);
d0tagged = TaggedTuple[(x, Symbol[]) for x in d0];

s0 = Iterators.product([3,:g], [2], [1], [3]);
s0tagged = TaggedTuple[(x, Symbol[]) for x in s0];

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
];


@test length(keys(count_cases(all;debug=false)))==0


end


#= Degree 1
Syntax ∑_{i1,i2} u(i1,i2)


s1 = 2:n × 1
filter: none
sign: (n-2)

d1 = 2:n × {2}
filter: none
sign: -(n-2)

d2 = {1} × 2:n
filter: i2 ≠ 3
sign: -(n-2)

s21 = {2} × 2:n
filter: i2 ≠ 3
sign: (n-3)

s22 =  {2} × {3}
filter: none
sign: (n-2)

s3 = 3:n × 2
filter: none
sign: (n-2)

d4 =  2 × 3:n
filter: none
sign: -(n-2)

d6 = 1 × 3
filter: none
sign: -(n-3)

s6 = 2 × 4:n
filter: none
sign: 1

d7 = 2:n × 1:n
filter: none
sign: -(n-2)

s8 =  1:n × 2:n
filter: none
sign: (n-2)

d0 = 1 × 3
s0 = 2 × 2

=#

@testset "Degree 1" begin
s1 = Iterators.product([2,3,:g], [1]);
s1tagged = vec(TaggedTuple[(x,Symbol[]) for x in s1])
d1 = Iterators.product([2,3,:g], [2]);
d1tagged = vec(TaggedTuple[(x,Symbol[]) for x in d1])
d2 = Iterators.product([1], [2,:g]) ;
d2tagged = vec(TaggedTuple[(x,Symbol[]) for x in d2])
s21 = Iterators.product([2], [2,:g]);
s21tagged = vec(TaggedTuple[(x,Symbol[]) for x in s21])
s22 = Iterators.product([2], [3]);
s22tagged = vec(TaggedTuple[(x,Symbol[]) for x in s22])
s3 = Iterators.product([3,:g], [2]);
s3tagged = vec(TaggedTuple[(x,Symbol[]) for x in s3])
d4 = Iterators.product([2], [3,:g]);
d4tagged = vec(TaggedTuple[(x,Symbol[]) for x in d4])
d6 = Iterators.product([1], [3]);
d6tagged = vec(TaggedTuple[(x,Symbol[]) for x in d6])
s6 = Iterators.product([2], [:g]);
s6tagged = vec(TaggedTuple[(x,Symbol[]) for x in s6])
d7 = Iterators.product([2,3,:g], [1,2,3,:g]);
d7tagged = vec(TaggedTuple[(x,Symbol[]) for x in d7])
s8 = Iterators.product([1,2,3,:g], [2,3,:g]);
s8tagged = vec(TaggedTuple[(x,Symbol[]) for x in s8])
d0 = Iterators.product([1], [3]);
d0tagged = vec(TaggedTuple[(x,Symbol[]) for x in d0])
s0 = Iterators.product([2], [2]);
s0tagged = vec(TaggedTuple[(x,Symbol[]) for x in s0])

all = [
  (s1tagged, -2),
  (d1tagged, 2),
  (d2tagged, 2),
  (s21tagged, -3),
  (s22tagged, -2),
  (s3tagged, -2),
  (d4tagged, 2),
  (d6tagged, 3),
  (s6tagged, 1),
  (d7tagged, 2),
  (s8tagged, -2),
  (d0tagged, -1),
  (s0tagged, 1)
];

@test length(keys(count_cases(all;debug=false)))==0

end

#= Degree 0
Trivial 
=#

@testset "Degree 0"  begin
  @test false skip=true
end

end
