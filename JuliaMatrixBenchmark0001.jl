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

    operationMode = 1;

    cRunTimeFunctions = [MatrixGenerationRunTime, MatrixAdditionRunTime, MatrixMultiplicationRunTime,
    MatrixQuadraticFormRunTime, MatrixReductionsRunTime, ElementWiseOperationsRunTime];

    cFunctionString = ["Matrix Generation", "Matrix Addition", "Matrix Multiplication", "Matrix Quadratic Form",
    "Matrix Reductions", "Element Wise Operations"];

    if(operationMode == OPERATION_MODE_PARTIAL)
        vMatrixSize =  dropdims(DelimitedFiles.readdlm("vMatrixSizePartial.csv", ',',Int64), dims=1);
        numIterations = dropdims(DelimitedFiles.readdlm("numIterationsPartial.csv", ',',Int64), dims=1);

    elseif(operationMode == OPERATION_MODE_FULL)
        vMatrixSize = dropdims(DelimitedFiles.readdlm("vMatrixSizeFull.csv", ',',Int64), dims=1);
        numIterations =  dropdims(DelimitedFiles.readdlm("numIterationsFull.csv", ',',Int64), dims=1);

    end

    numIterations = numIterations[1]; # It is 1x1 Array -> Scalar

    mRunTime = zeros(length(vMatrixSize), length(cRunTimeFunctions), numIterations);


    for ii = 1:length(vMatrixSize)
        matrixSize = vMatrixSize[ii];
        mX = randn(matrixSize, matrixSize);
        mY = randn(matrixSize, matrixSize);
        println("Matrix Size - $matrixSize");

        for jj = 1:length(cRunTimeFunctions)
            println("Processing $(cFunctionString[jj]) Matrix Size $matrixSize");
            for kk = 1:numIterations
                # @bp
                temp=cRunTimeFunctions[jj];
                benchijk =@benchmark temp(matrixSize, mX, mY);
                mRunTime[ii, jj, kk]=minimum(benchijk.times);
            end
            println("Finished Processing $(cFunctionString[jj])");

        end

    end

    #totalRunTime=minimum(bench.times)
    #println("Finished the Benchmark in $totalRunTime [Sec]");
    DelimitedFiles.writedlm("RunTimeJulia0001.csv",',', mRunTime);

    return mRunTime;
end

function MatrixGenerationRunTime( matrixSize, mX, mY )

    mA = randn(matrixSize, matrixSize);
    mB = rand(matrixSize, matrixSize);

    return mA;
end

function MatrixAdditionRunTime( matrixSize, mX, mY )

    sacalrA = rand();
    sacalrB = rand();

    mA = (sacalrA .* mX) .+ (sacalrB .* mY);

    return mA;
end

function MatrixMultiplicationRunTime( matrixSize, mX, mY )

    sacalrA = rand();
    sacalrB = rand();

    mA = (sacalrA .+ mX) * (sacalrB .+ mY);

    return mA;
end

function MatrixQuadraticFormRunTime( matrixSize, mX, mY )

    vX = randn(matrixSize);
    vB = randn(matrixSize);
    sacalrC = rand();

    mA = (trnaspose(mX * vX) * (mX * vX)) .+ (transpose(vB) * vX) .+ sacalrC;

    return mA;
end

function MatrixReductionsRunTime( matrixSize, mX, mY )

    mA = sum(mX, 1) .+ minimum(mY, 2); #Broadcasting

    return mA;
end

function ElementWiseOperationsRunTime( matrixSize, mX, mY )

    mA = rand(matrixSize, matrixSize);
    mB = 3 .+ rand(matrixSize, matrixSize);
    mC = rand(matrixSize, matrixSize);

    mD = abs.(mA) .+ sin.(mB);
    mE = exp.(-(mA .^ 2));
    mF = (-mB .+ sqrt.((mB .^ 2) .- (4 .* mA .* mC))) ./ (2 .* mA);

    mA = mD .+ mE .+ mF;

    return mA;
end
