
clearconsole();

# ]add BenchmarkTools
# ]add MAT
using BenchmarkTools # for benchmark
BenchmarkTools.DEFAULT_PARAMETERS.samples = 700;
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1;
using DelimitedFiles # for readdlm, writedlm
using Statistics # for median
using LinearAlgebra # for eigen, pinv, svd, cholesky
using Random # for randperm
using MAT


operationMode = 2; # 0 for test only # 1 for partial benchmark # 2 for full benchmark


include("JuliaBench.jl");
tRunTime, mRunTime= JuliaBench(operationMode);

writedlm("RunTimeData\\RunTimeJulia$(BLAS.vendor())Table.csv", tRunTime,',');

file = matopen("RunTimeData\\RunTimeJulia$(BLAS.vendor()).mat", "w")
write(file, "mRunTime", mRunTime)
close(file)


include("JuliaBenchSIMD.jl")
tRunTime, mRunTime= JuliaBenchSIMD(operationMode);

writedlm("RunTimeData\\RunTimeJulia$(BLAS.vendor())SIMDTable.csv", tRunTime,',');

file = matopen("RunTimeData\\RunTimeJulia$(BLAS.vendor())SIMD.mat", "w")
write(file, "mRunTime", mRunTime)
close(file)
