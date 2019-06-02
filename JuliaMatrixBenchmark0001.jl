# ----------------------------------------------------------------------------------------------- #
# Julia Matrix Operations Benchmark - Test Suite 0001
# Reference:
#   1. C.
# Remarks:
#   1.  W.
# TODO:
#   1.  A
#   Release Notes:
#   -   1.0.002     10/02/2017  Royi Avital
#       *   Added generation of 'mX' and 'mY' once outside the functions.
#       *   Fixed issue with the Quadratic Form.
#       *   Optimized creation of scalars and vectors.
#   -   1.0.001     09/02/2017  Royi Avital
#       *   Added 'MatrixExpRunTime()' and 'MatrixSqrtRunTime()'.
#       *   Added Quadratic Matrix Form Calculation 'MatrixQuadraticFormRunTime()'.
#       *   Added Univariate Quadratic Function Root to 'ElementWiseOperationsRunTime()'.
#       *   Updated 'MatrixGenerationRunTime()' to include Uniform Random Number Generation.
#       *   Fixed issue with 'CalcDistanceMatrixRunTime'.
#   -   1.0.000     09/02/2017  Royi Avital
#       *   First release version.
# ----------------------------------------------------------------------------------------------- #

function JuliaMatrixBenchmark0001( operationMode = 2 )


    OPERATION_MODE_PARTIAL  = 1; # For Testing (Runs Fast)
    OPERATION_MODE_FULL     = 2;

    cRunTimeFunctions = [MatrixGenerationRunTime, MatrixAdditionRunTime, MatrixMultiplicationRunTime,
    MatrixQuadraticFormRunTime, MatrixReductionsRunTime, ElementWiseOperationsRunTime];

    cFunctionString = ["Matrix Generation", "Matrix Addition", "Matrix Multiplication", "Matrix Quadratic Form",
    "Matrix Reductions", "Element Wise Operations"];

    if(operationMode == OPERATION_MODE_PARTIAL)
        vMatrixSize = round.(Int64, dropdims(DelimitedFiles.readdlm("vMatrixSizePartial.csv", ','), dims=1));
        numIterations = round.(Int64, dropdims(DelimitedFiles.readdlm("numIterationsPartial.csv", ','), dims=1));
    elseif(operationMode == OPERATION_MODE_FULL)
        vMatrixSize = round.(Int64, dropdims(DelimitedFiles.readdlm("vMatrixSizeFull.csv", ','), dims=1));
        numIterations = round.(Int64, dropdims(DelimitedFiles.readdlm("numIterationsFull.csv", ','), dims=1));
    end

    numIterations = numIterations[1]; # It is 1x1 Array -> Scalar

    mRunTime = zeros(length(vMatrixSize), length(cRunTimeFunctions), numIterations);

    #bench=@benchmark begin
    for ii = 1:length(vMatrixSize)
        matrixSize = vMatrixSize[ii];
        mX = randn(matrixSize, matrixSize);
        mY = randn(matrixSize, matrixSize);
        println("Matrix Size - $matrixSize");
        for jj = 1:length(cRunTimeFunctions)
            println("Processing $(cFunctionString[jj]) Matrix Size $matrixSize");
            for kk = 1:numIterations
                mA, mRunTime[ii, jj, kk] = cRunTimeFunctions[jj](matrixSize, mX, mY);
            end
            println("Finished Processing $(cFunctionString[jj])");
        end
    end
    #end
    #totalRunTime=minimum(bench.times)

    #mRunTime = median(bench.times, 3);
    #mRunTime = squeeze(mRunTime, 3);

    #println("Finished the Benchmark in $totalRunTime [Sec]");

    #writecsv("RunTimeJulia0001.csv", mRunTime);
    return 0;

end

function MatrixGenerationRunTime( matrixSize, mX, mY )

    bench=@benchmark begin
        mA = randn(matrixSize, matrixSize);
        mB = rand(matrixSize, matrixSize);
    end
    runTime=minimum(bench.times)
    return mA, runTime;
end

function MatrixAdditionRunTime( matrixSize, mX, mY )

    sacalrA = rand();
    sacalrB = rand();

    bench=@benchmark begin
        mA = (sacalrA .* mX) .+ (sacalrB .* mY);
    end
    runTime=minimum(bench.times)

    return mA, runTime;
end

function MatrixMultiplicationRunTime( matrixSize, mX, mY )

    sacalrA = rand();
    sacalrB = rand();

    bench=@benchmark begin
        mA = (sacalrA .+ mX) * (sacalrB .+ mY);
    end
    runTime=minimum(bench.times)

    return mA, runTime;
end

function MatrixQuadraticFormRunTime( matrixSize, mX, mY )

    vX = randn(matrixSize);
    vB = randn(matrixSize);
    sacalrC = rand();

    bench=@benchmark begin
        mA = (trnaspose(mX * vX) * (mX * vX)) .+ (transpose(vB) * vX) .+ sacalrC;
    end
    runTime=minimum(bench.times)

    return mA, runTime;
end

function MatrixReductionsRunTime( matrixSize, mX, mY )

    bench=@benchmark begin
        mA = sum(mX, 1) .+ minimum(mY, 2); #Broadcasting
    end
    runTime=minimum(bench.times)

    return mA, runTime;
end

function ElementWiseOperationsRunTime( matrixSize, mX, mY )

    mA = rand(matrixSize, matrixSize);
    mB = 3 .+ rand(matrixSize, matrixSize);
    mC = rand(matrixSize, matrixSize);

    bench=@benchmark begin
        mD = abs.(mA) .+ sin.(mB);
        mE = exp.(-(mA .^ 2));
        mF = (-mB .+ sqrt.((mB .^ 2) .- (4 .* mA .* mC))) ./ (2 .* mA);
    end
    runTime=minimum(bench.times)

    mA = mD .+ mE .+ mF;

    return mA, runTime;
end
