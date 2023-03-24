using Plots
using LaTeXStrings
using SpecialFunctions
using CSV, DataFrames
using QuadGK

function test_func(x, y)
    x*y*abs(sin(x*y^3))*exp(-x^2/y)
end

length_xval = 100
length_yval = 100

results_file = open("Results.csv", "w")
IOStream("Results.csv")
write(results_file, "x, y, x/y, f\n")

Threads.@threads for (i,j) in collect(Iterators.product(range(1,10,length_xval), range(1,10,length_yval)))
    ratio = i/j
    value = test_func(i, j)
    write(results_file, join((i,j,ratio, value),","),"\n")
end

close(results_file)