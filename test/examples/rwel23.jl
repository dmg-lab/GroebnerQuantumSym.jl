function rwel23_data(n::Int)
#Put requirements for n here
@assert n > 3

G1 = g1_named(n)
u = magic_unitary(n)
E1 = extra_relations(u)
QuantumGB.add!(E1, G1)


to_be_reduced()=rwel(2,3; u=u)
free_vars= ()



reduction() = begin (-sum([E1["rinjcs$i"] for i in 2:n])
  +sum([E1[Symbol("rwelcs$i")] for i in 2:n if i != 3])
  +sum([E1["rcs$(i)2"] for i in 3:n])
  -sum([E1["rrs2$(j)"] for j in 3:n])
  -sum([(E1["rrs$i$j"]-E1["rcs$i$j"]) for i in 3:n for j in 3:n])
  +sum([E1["rwel$(i)3"] for i in 4:n])
  +(n-2)*sum([E1["rs$i"] for i in 2:n])
  -(n-2)*sum([E1["cs$i"] for i in 2:n])
)
end

  return to_be_reduced, free_vars, reduction
end

@testset "rwel23" verbose=true begin
  println("Testing examples of rwel23...")
  @testset "n=4:15" begin
    for n in 4:15
      @testset "n=$n" begin
        to_be_reduced, free_vars, reduction = rwel23_data(n)
        @test to_be_reduced()+reduction() == 0
      end
    end
  end
  @testset "n=35" begin
    to_be_reduced, free_vars, reduction = rwel23_data(20)
    @test to_be_reduced()+reduction() == 0
  end
end


