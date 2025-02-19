using QuantumGB
using ProgressMeter
#=
include("../examples/bg223t_data.jl")

x = 7
ns = collect(6:10)
@test_reduction ns "../examples/bg223t_data.jl"
@test_reduction x "../examples/bg223t_data.jl"

@test_reduction x reduction to_be_reduced free_vars
@test_reduction ns reduction to_be_reduced free_vars
=#

macro _test_core(reduction::Symbol,to_be_reduced::Symbol,free_vars::Symbol)
  println("ello")
  println(eval(reduction))
  vars = Meta.parse.((map(x->String(x)[1:end-1],keys(eval(free_vars)))))
  
  free_vars_expr = eval(free_vars)
  loops = Expr(:block)
  reduction_expr = Meta.parse(eval(reduction))
  to_be_reduced_expr = Meta.parse(eval(to_be_reduced))

  loops = quote
    local _reduction = $reduction_expr
    local _to_be_reduced = $to_be_reduced_expr
    local red = _to_be_reduced+_reduction
    iszero(red) && continue
    push!(report,(reduction_string(G1,red),red))
  end

  for (key, value) in zip(vars,Meta.parse.(values(free_vars_expr)))
    loop = :(for $key in $value
      $(loops.args...)
    end)
    loops = Expr(:block, loop)
  end

  quote
    local report = []

    $(loops.args...)

    report
  end
end

macro test_reduction(n::Symbol)
  _n_expr = eval(n)

  if typeof(_n_expr) == Int64
    expr = quote
      local n = $n 
      local G1 = g1_named(n)
      local u = magic_unitary(G1)
      local reduction, to_be_reduced, free_vars = data_f(n)

      report = @_test_core($reduction,to_be_reduced,free_vars)

      length(report) == 0 
    end
  else
    expr = quote
      local rep = true

      for n in $_n_expr
        local G1 = g1_named(n)
        local u = magic_unitary(G1)
        local reduction, to_be_reduced, free_vars = $data_f(n)
        println(reduction)

        report = @_test_core(reduction,to_be_reduced,free_vars)

        length(report) == 0 && continue 
        println("Failed for n = ",n)
        println("Reduction: ")
        println.(report) 
        rep = false
        end
        return rep
      end
    end
    return expr
end

macro test_reduction(ns::Symbol,path::String)
  include(path)
  quote 
    @test_reduction ns
  end
end
