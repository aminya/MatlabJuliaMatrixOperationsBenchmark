function AnalysisJuliaPlotter()
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
    ii=1;
    for fun in sFunNameMatlab

        plt=plot( vMatrixSizeMatlab,[ mRunTimeMatlab[ii,:], mRunTimeJuliamkl[ii,:] ],
        labels=["MATLAB" "Julia-MKL"],legend=:bottomright, markershape =:auto,markersize=2,
        xlabel="Matrix Size", ylabel="Run Time  [micro Seconds]",
        title="$fun", xscale=:log10, yscale=:log10,dpi=300 );

        plotJuliaSIMD=occursin.(fun,sFunNameJuliamklSIMD); # if 1 will plot JuliaSIMD
        if any(plotJuliaSIMD)
            plt=plot!(vMatrixSizeJuliamklSIMD,dropdims(mRunTimeJuliamklSIMD[plotJuliaSIMD,:],dims=1),
            label="Julia-MKL-SIMD",markershape =:auto,markersize=2);
        end
        display(plt); # display in plot pane
        if(generateImages == 1)
            savefig("Figures\\Julia\\Figure$(ii).png");
        end

        ii=ii+1;
    end

end