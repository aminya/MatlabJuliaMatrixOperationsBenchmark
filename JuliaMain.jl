
clearconsole();

using BenchmarkTools # for benchmark
using DelimitedFiles # for readdlm, writedlm
using Statistics # for median
using LinearAlgebra # for eigen, pinv, svd, cholesky
using Random # for randperm

include("JuliaBench.jl");
# include("JuliaBenchOptimized.jl")

operationMode = 2; # 0 for test only # 1 for partial benchmark # 2 for full benchmark 

tRunTime= JuliaBench(operationMode);
# tRunTime=JuliaBenchOptimized(operationMode);

# for debugging:
 # using JuliaInterpreter
 # @interpret JuliaBench(operationMode)
# use @bp in different parts of the code for breakpoint
