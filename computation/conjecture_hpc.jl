#Code to Run on the cluster
using QuantumGB

#n will be set by *.job file

deg_bound = 8
str = check_conjecture(n; deg_bound=deg_bound)

# Write the output to a file
filepath = "../data/conjecture_n_$(n).txt"

open(filepath, "w") do io
    println(io, str)
end


