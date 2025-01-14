export ColoredBlock, ColoredPartition

export colored_block, cb,  dim, max_value, cube
export colored_partition

struct ColoredBlock
  block::Vector{Int}
  color::Int
end

function colored_block(b::Vector{Int}, color::Int)
  return ColoredBlock(b, color)
end

cb(b::Vector{Int}, color::Int) = colored_block(b, color)
cb(b::Int, color::Int) = colored_block([b], color)

Base.:(==)(cb1::ColoredBlock, cb2::ColoredBlock) = cb1.block == cb2.block
Base.isequal(cb1::ColoredBlock, cb2::ColoredBlock) = cb1 == cb2
Base.hash(cb::ColoredBlock, h::UInt) = hash(cb.block, hash(cb.color, h))
Base.isless(cb1::ColoredBlock, cb2::ColoredBlock) = cb1.block < cb2.block
function Base.:(<=)(cb1::ColoredBlock, cb2::ColoredBlock) 
  @assert dim(cb1) == dim(cb2)
  return all([cb1.block[i] <= cb2.block[i] for i in 1:dim(cb1)])
end

function Base.:(<)(cb1::ColoredBlock, cb2::ColoredBlock)
  @assert dim(cb1) == dim(cb2)
  return all([cb1.block[i] < cb2.block[i] for i in 1:dim(cb1)])
end

Base.show(io::IO, cb::ColoredBlock) =  print(io, cb.block, _index_number(cb.color))

function Base.:(+)(a::ColoredBlock, b::ColoredBlock)
  @assert a.block == b.block
  return colored_block(a.block, a.color + b.color)
end

function add_color(cb::ColoredBlock, color::Int)
  return colored_block(cb.block, cb.color + color)
end

dim(cb::ColoredBlock) = length(cb.block)
max_value(cb::ColoredBlock) = maximum(cb.block)
min_value(cb::ColoredBlock) = minimum(cb.block)

function cube(cb::ColoredBlock) 
  S = cb.block
  d = dim(cb)
  # Create a vecor of length S, with each possible combination of values [0,1}
  corners = [bitstring(i)[end-d+1:end] for i in 1:(2^d-1)]
  corners = [-1 * parse.(Int, collect(corner)) + S for corner in corners]
  return corners 
end

#=
--[1]-- should be the block [1,..n)

x = cb(1, 1)
y = cb(2, 2)
z = cb(3, 1)

z2 = cb([2,2,2], 1)
cube(z2)

x == y

Y = cb([2,1], 1)
X = cb([1,2], 1)
X >= Y


x * y

Create a d dimensional array of size n^d
n = 3  # size of each dim
d = 4  # number of dimensions
X = Array{UInt, 4}(undef, ntuple(_ -> 3, 4)...)
X[1,1,5,1] = 1

=#

mutable struct ColoredPartition
  block::ColoredBlock
  children::Vector{Union{ColoredPartition, Nothing}}
  parents::Vector{Union{ColoredPartition, Nothing}}
  function ColoredPartition(block::ColoredBlock, children::Vector{Union{ColoredPartition, Nothing}}, parents::Vector{Union{ColoredPartition, Nothing}})
    @assert all([dim(block) == dim(child.block) for child in children if child isa ColoredPartition])
    @assert length(parents) == 2
    return new(block, children, parents)
  end
end

function colored_partition(block::ColoredBlock)
  return ColoredPartition(block, Vector{Union{ColoredPartition, Nothing}}(), Union{ColoredPartition, Nothing}[nothing,nothing])
end

function add_child!(cp::ColoredPartition, child::ColoredPartition)
  push!(cp.children, child)
  #replace first nothing with parent in child.parents
  if isnothing(child.parents[1])
    child.parents[1] = cp
  elseif isnothing(child.parents[2])
    child.parents[2] = cp
  else
    error("Child already has two parents")
  end
end


#=
colored_partition(cb(1, 1))

=#
