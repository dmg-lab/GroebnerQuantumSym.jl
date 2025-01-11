export ColoredBlock, ColoredPartition

export colored_block, cb, cpp, gc_refinement

mutable struct ColoredBlock
  block::Int
  color::Int
end

function colored_block(b::Int, color::Int)
  return ColoredBlock(b, color)
end

cb(b::Int, color::Int) = colored_block(b, color)


Base.:(==)(cb1::ColoredBlock, cb2::ColoredBlock) = cb1.block == cb2.block
Base.isequal(cb1::ColoredBlock, cb2::ColoredBlock) = cb1 == cb2
Base.hash(cb::ColoredBlock, h::UInt) = hash(cb.block, hash(cb.color, h))
Base.isless(cb1::ColoredBlock, cb2::ColoredBlock) = cb1.block < cb2.block
Base.:(<=)(cb1::ColoredBlock, cb2::ColoredBlock) = cb1.block <= cb2.block

Base.show(io::IO, cb::ColoredBlock) =  print(io, cb.block, _index_number(cb.color))

function Base.:(+)(a::ColoredBlock, b::ColoredBlock)
  @assert a.block == b.block
  return colored_block(a.block, a.color + b.color)
end

function add_color(cb::ColoredBlock, color::Int)
  return colored_block(cb.block, cb.color + color)
end





#=
--[1]-- should be the block [1,..n)

x = cb(1, 1)
y = cb(2, 2)
z = cb(3, 1)

x == y

S1 = Set([y, x])
S2 = Set([z, y])

x * y

=#

mutable struct ColoredPartition
  blocks::Vector{ColoredBlock}
  function ColoredPartition(S::Vector{ColoredBlock})

  @assert unique(S) == S
    @assert length(S) > 0

    new_block = Vector{ColoredBlock}()

    if !issorted(S)
      S = sort(S)
    end
    # handle the coarsening in the constructor
    for i in 1:length(S)
      if i == 1 
        push!(new_block, S[i])
        continue
      end
      S[i].color == S[i-1].color && continue
      push!(new_block, S[i])
    end

    return new(new_block)
  end
end

#=
Partition 1: [1,3), [4,7), [7,n]
Partition 2: [1,2), [2,5), [5,n]
CP1 = cpp([cb(1,1), cb(4,2), cb(7,2)])
CP2 = cpp([cb(1,1), cb(2,2), cb(5,1)])
gc_refinement(CP1, CP2)
Answer should be [1,2)(2), [2,4)(3), [4,5)(4), [5,n)(3)

get_next_smallest_block(CP2, cb([4],2))
CP1 + CP2
CP1 == CP1
CP2
dir_product(CP1,CP2)
=#

function Base.:(==)(cp1::ColoredPartition, cp2::ColoredPartition)
  length(cp1.blocks) == length(cp2.blocks) || return false
  for (b1, b2) in zip(cp1.blocks, cp2.blocks)
    if b1 != b2
      return false
    end
  end
  return true
end

function Base.hash(cp::ColoredPartition, h::UInt)
  hash(cp.blocks, h)
end


function Base.show(io::IO, cp::ColoredPartition) 
  S1 = sort(collect(cp.blocks))
  ans = ""
  for cb in S1
    ans *= string(cb) * " "
  end
  print(io, ans)
end

colored_partition(b::Set{ColoredBlock}) = ColoredPartition(b)
cpp(b::Vector{ColoredBlock}) = ColoredPartition(b)

function get_next_smallest_block(cp::ColoredPartition, b::ColoredBlock)
  curr = nothing
  for cb in cp.blocks
    if cb <= b
      curr = cb
    else
      return curr
    end
  end
  return curr
end

function gc_refinement(cp1::ColoredPartition, cp2::ColoredPartition)
  blocks = Vector{ColoredBlock}()
  function _refine(cp1,cp2)
  S1 = cp1.blocks
  S2 = cp2.blocks
    for b1 in S1
      b2 = get_next_smallest_block(cp2, b1)
      if isnothing(b2) 
        push!(blocks, b1)
      elseif b1 == b2
        push!(blocks, cb(b1.block, b1.color + b2.color))
        filter!(x -> x != b2, S2)
      else
        push!(blocks, cb(b1.block, b1.color + b2.color))
      end
    end
  end
  _refine(cp1, cp2) 
  _refine(cp2, cp1)
  return ColoredPartition(sort(blocks))
end


function dir_product(cp1::ColoredPartition, cp2::ColoredPartition)
  blocks = Set{ColoredBlock}()
  S1 = sort(collect(cp1.blocks))
  S2 = sort(collect(cp2.blocks))
  for b1 in S1
    for b2 in S2
      push!(blocks, b1 * b2)
    end
  end
  return ColoredPartition(blocks)
end


Base.:(+)(cp1::ColoredPartition, cp2::ColoredPartition) = gc_refinement(cp1, cp2)


mutable struct ComplexColoredBlock
  block::Int
  color::Int
  partition::ColoredPartition
end

function complex_colored_block(b::Int, color::Int, partition::ColoredPartition)
  return ComplexColoredBlock(b, color, partition)
end

ccb(b::Int, color::Int, partition::ColoredPartition) = complex_colored_block(b, color, partition)


function Base.:(==)(cb1::ComplexColoredBlock, cb2::ComplexColoredBlock)
    for (b1, b2) in zip(cb1.block.blocks, cb2.block.blocks)
      if b1 != b2
        return false
      end
    end 
end
Base.isequal(cb1::ColoredBlock, cb2::ColoredBlock) = cb1 == cb2
Base.hash(cb::ColoredBlock, h::UInt) = hash(cb.block, hash(cb.color, h))
Base.isless(cb1::ColoredBlock, cb2::ColoredBlock) = cb1.block < cb2.block
Base.:(<=)(cb1::ColoredBlock, cb2::ColoredBlock) = cb1.block <= cb2.block

Base.show(io::IO, cb::ColoredBlock) =  print(io, cb.block, _index_number(cb.color))

function Base.:(+)(a::ColoredBlock, b::ColoredBlock)
  @assert a.block == b.block
  return colored_block(a.block, a.color + b.color)
end

function add_color(cb::ColoredBlock, color::Int)
  return colored_block(cb.block, cb.color + color)
end

