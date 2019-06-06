
clearconsole();

using DelimitedFiles
using BenchmarkTools
using LinearAlgebra
using Random

include("JuliaBench.jl");
# include("JuliaBenchOptimized.jl")

operationMode = 1; # 1 for partial fast benchmark # 2 for full benchmark

mRunTime= JuliaBench(operationMode);
# mRuntime=JuliaBenchOptimized(operationMode);

# for debugging:
# using JuliaInterpreter
# @interpret JuliaFuns(operationMode)
# use @bp in different parts of the code for breakpoint
