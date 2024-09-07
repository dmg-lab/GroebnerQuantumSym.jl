#Code to Run on the cluster

#n will be set by *.job file

str = check_conjecture(n)

# Write the output to a file
filepath = "../data/conjecture_n_$(n).txt"

open(filepath, "w") do io
    println(io, str)
end


