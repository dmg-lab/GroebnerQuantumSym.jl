using Test
using QuantumGB

println("Testing...")


@testset "QuantumGB" verbose=true begin


@testset "Generic" verbose=true begin
for file in readdir(joinpath(@__DIR__, "..", "test", "generic"))
    if endswith(file, ".jl")
        include(joinpath(@__DIR__, "..", "test", "generic", file))
    end
end
end

@testset "Examples" verbose=true begin
for file in readdir(joinpath(@__DIR__, "..", "test", "examples"))
  if endswith(file, ".jl")
    include(joinpath(@__DIR__, "..", "test", "examples", file))
  end
end
end

end;

