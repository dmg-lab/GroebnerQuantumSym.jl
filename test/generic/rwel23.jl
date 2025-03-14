
@testset "rwel23" begin
println("""
Testing the general case of rwel23...""")


@testset "Degree 2" begin

L1 = lexicon([[1],[2],[3],[:g]])
L2 = lexicon([[1],[2],[3],[:g]])
L3 = lexicon([[1],[2],[3],[:g],[Symbol("g⩓≠i1")]])
L4  = lexicon([[1],[2],[3],[:g],[Symbol("g⩓≠i2")]])
L = L1 * L2 * L3 * L4

s1 = create_preimage_in(L, [[2,3,:g],[2],[2,3,:g],[3,:g]], 1; filter=[(3,1)]);
d1 = create_preimage_in(L, [[2,3,:g],[3,:g],[2,3,:g],[1]], -1; filter=[(3,1)]);
d2 = create_preimage_in(L, [[2],[2,3,:g],[3,:g],[2,:g]], -1; filter=[(4,2)]);
s2 = create_preimage_in(L, [[3,:g],[2,3,:g],[1],[2,:g]], 1; filter=[(4,2)]);
d3 = create_preimage_in(L, [[3,:g],[2],[1,2,3,:g],[2,3,:g]], -1; filter=[(3,1),(4,2)]);
s4 = create_preimage_in(L, [[2],[3,:g],[2,3,:g],[1,2,3,:g]], 1; filter=[(3,1),(4,2)]);
s5 = create_preimage_in(L, [[3,:g],[3,:g],[2,3,:g],[1,2,3,:g]], 1; filter=[(3,1),(4,2)]);
d6 = create_preimage_in(L, [[3,:g],[3,:g],[1,2,3,:g],[2,3,:g]], -1; filter=[(3,1),(4,2)]);
d7 = create_preimage_in(L, [[2],[:g],[3,:g],[3]], -1);
s7 = create_preimage_in(L, [[3,:g],[:g],[1],[3]], 1);
d0 = create_preimage_in(L, [[2],[2],[3,:g],[3]], -1);
s0 = create_preimage_in(L, [[3,:g],[2],[1],[3]], 1);

S = s1 + s2 + s4 + s5 + s7 + s0 + d1 + d2 + d3 + d6 + d7 + d0
@test iszero(S)
end
@testset "Degree 1" begin

L1 = lexicon([[1],[2],[3],[:g]])
L2 = lexicon([[1],[2],[3],[:g]])
L = L1 * L2
R = L.coefficient_ring
n = gen(R)

s1 = create_preimage_in(L, [[2,3,:g],[1]], (n-2));
d1 = create_preimage_in(L, [[2,3,:g],[2]], -(n-2));
d2 = create_preimage_in(L, [[1],[2,:g]], -(n-2));
s21 = create_preimage_in(L, [[2],[2,:g]], (n-3));
s22 = create_preimage_in(L, Vector{Any}[[2],[3]], (n-2));
s3 = create_preimage_in(L, [[3,:g],[2]], (n-2));
d4 = create_preimage_in(L, [[2],[3,:g]], -(n-2));
d6 = create_preimage_in(L, Vector{Any}[[1],[3]], -(n-3));
s6 = create_preimage_in(L, [[2],[:g]], 1);
d7 = create_preimage_in(L, [[2,3,:g],[1,2,3,:g]], -(n-2));
s8 = create_preimage_in(L, [[1,2,3,:g],[2,3,:g]], (n-2));
d0 = create_preimage_in(L, Vector{Any}[[1],[3]], -1);
s0 = create_preimage_in(L, Vector{Any}[[2],[2]], 1);

S = s1 + s21 + s22 + s3 + s6 + s8 + s0 + d1 + d2 + d4 + d6 + d7 + d0
@test iszero(S)

end

@testset "Degree 0"  begin
  @test false skip=true
end
end


