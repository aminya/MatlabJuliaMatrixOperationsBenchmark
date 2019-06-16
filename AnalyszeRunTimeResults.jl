# using Pkg
# Pkg.add("Plots")
# Pkg.add("GR")

using DelimitedFiles
using LinearAlgebra
using Plots
# Setting Enviorment Parameters
generateImages=1;

# Loading Data
tRunTimeMatlab = readdlm("RunTimeData\\RunTimeMatlabTable.csv", ',');
mRunTimeMatlab=tRunTimeMatlab[2:end,2:end];
vMatrixSizeMatlab=tRunTimeMatlab[1,2:end];
sFunNameMatlab=tRunTimeMatlab[2:end,1];

tRunTimeJuliamkl = readdlm("RunTimeData\\RunTimeJuliamklTable.csv", ',');
mRunTimeJuliamkl=tRunTimeJuliamkl[2:end,2:end];
vMatrixSizeJuliamkl=tRunTimeJuliamkl[1,2:end];
sFunNameJuliamkl=tRunTimeJuliamkl[2:end,1];

tRunTimeJuliamklSIMD = readdlm("RunTimeData\\RunTimeJuliamklSIMDTable.csv", ',');
mRunTimeJuliamklSIMD=tRunTimeJuliamklSIMD[2:end,2:end];
vMatrixSizeJuliamklSIMD=tRunTimeJuliamklSIMD[1,2:end];
sFunNameJuliamklSIMD=tRunTimeJuliamklSIMD[2:end,1];

# Displaying Results
for ii = 1:size(sFunNameMatlab,2)

    plot( vMatrixSizeMatlab,[ mRunTimeMatlab[ii,:], mRunTimeJuliamkl[ii,:] ],
    labels=["MATLAB" "Julia-MKL"],
    xlabel="Matrix Size", ylabel="Run Time  [micro Seconds]",
    title="$(sFunNameMatlab[ii])", xaxis=:log, yaxis=:log );

    plotJuliaSIMD=ismember( fun, sFunNameJuliamklSIMD ); # if 1 will plot JuliaSIMD
    if any(plotJuliaSIMD)
        plot!(vMatrixSizeJuliamklSIMD,mRunTimeJuliamklSIMD[plotJuliaSIMD,:],
        label="Julia-MKL-SIMD",
        xaxis=:log, yaxis=:log,);
    end

    if(generateImages == 1)
        savefig("Figure$ii.png");
    end

end
