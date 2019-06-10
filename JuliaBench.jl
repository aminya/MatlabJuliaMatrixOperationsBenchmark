
function JuliaBench(operationMode)

    allFunctions = [MatrixGeneration, MatrixAddition, MatrixMultiplication, MatrixQuadraticForm, MatrixReductions, ElementWiseOperations, MatrixExp, MatrixSqrt, Svd, Eig, CholDec, MatInv, LinearSystem, LeastSquares, CalcDistanceMatrix, KMeans];


    allFunctionsString = ["Matrix Generation", "Matrix Addition", "Matrix Multiplication", "Matrix Quadratic Form", "Matrix Reductions", "Element Wise Operations", "Matrix Exponential", "Matrix Square Root", "Singular Value Decomposition", "Eigen Decomposition","Cholesky Decomposition", "Matrix Inversion", "Linear System Solution", "Linear Least Squares", "Squared Distance Matrix", "K-Means"];

    if(operationMode == 1) # partial fast benchmark
        vMatrixSize =  dropdims(DelimitedFiles.readdlm("vMatrixSizePartial.csv", ',',Int64), dims=1);
        numIterations = dropdims(DelimitedFiles.readdlm("numIterationsPartial.csv", ',',Int64), dims=1);

    elseif(operationMode == 2) # full fast benchmark
        vMatrixSize = dropdims(readdlm("vMatrixSizeFull.csv", ',',Int64), dims=1);
        numIterations =  dropdims(readdlm("numIterationsFull.csv", ',',Int64), dims=1);

    end

    numIterations = numIterations[1]; # It is 1x1 Array -> Scalar

    mRunTime = zeros(length(vMatrixSize), length(allFunctions), numIterations);


    for ii = 1:length(vMatrixSize)
        matrixSize = vMatrixSize[ii];
        mX = randn(matrixSize, matrixSize);
        mY = randn(matrixSize, matrixSize);
        println("Matrix Size - $matrixSize");

        jj=1;
        for fun in allFunctions
            println("Processing $(allFunctionsString[jj]) - MatrixSize= $matrixSize");

            for kk = 1:numIterations
                benchmarkableIJK =@benchmarkable $fun($matrixSize, $mX, $mY);
                tune!(benchmarkableIJK)
                benchmarkableIJK
                # benchIJK =@benchmark $fun($matrixSize, $mX, $mY);
                # mRunTime[ii, jj, kk]=minimum(benchijk.times);

            end
            jj+=1;
        end

    end

    writedlm("RunTimeJulia0001.csv", mRunTime);

    return mRunTime;
end

function MatrixGeneration( matrixSize, mX, mY )

    mA = randn(matrixSize, matrixSize);
    mB = rand(matrixSize, matrixSize);

    return mA;
end

function MatrixAddition( matrixSize, mX, mY )

    scalarA = rand();
    scalarB = rand();

    mA = (scalarA .* mX) .+ (scalarB .* mY);

    return mA;
end

function MatrixMultiplication( matrixSize, mX, mY )

    scalarA = rand();
    scalarB = rand();

    mA = (scalarA .+ mX) * (scalarB .+ mY);

    return mA;
end

function MatrixQuadraticForm( matrixSize, mX, mY )

    vX = randn(matrixSize);
    vB = randn(matrixSize);
    scalarC = rand();

    mA = (transpose(mX * vX) * (mX * vX)) .+ (transpose(vB) * vX) .+ scalarC;

    return mA;
end

function MatrixReductions( matrixSize, mX, mY )

    mA = sum(mX, dims=1) .+ minimum(mY, dims=2); #Broadcasting

    return mA;
end

function ElementWiseOperations( matrixSize, mX, mY )

    mA = rand(matrixSize, matrixSize);
    mB = 3 .+ rand(matrixSize, matrixSize);
    mC = rand(matrixSize, matrixSize);

    mD = abs.(mA) .+ sin.(mB);
    mE = exp.(-(mA .^ 2));
    mF = (-mB .+ sqrt.((mB .^ 2) .- (4 .* mA .* mC))) ./ (2 .* mA);

    mA = mD .+ mE .+ mF;

    return mA;
end

function MatrixExp( matrixSize, mX, mY )

    mA = exp(mX);

    return mA;
end

function MatrixSqrt( matrixSize, mX, mY  )

    mY = transpose(mX) * mX;

    mA = sqrt(mY);

    return mA;
end

function Svd( matrixSize, mX, mY  )

    F = svd(mX, full = false);

    return mA=0;
end

function Eig( matrixSize, mX, mY  )

    mE = eigvals(mX);

    return mA=0;
end

function CholDec( matrixSize, mX, mY  )

    mY = transpose(mX) * mX;

    mA = cholesky(mY);

    return mA;
end

function MatInv( matrixSize, mX, mY  )

    mY = transpose(mX) * mX;

    mA = inv(mY);
    mB = pinv(mX);

    mA = mA .+ mB;

    return mA;
end


function LinearSystem( matrixSize, mX, mY  )

    mB = randn(matrixSize, matrixSize);
    vB = randn(matrixSize);

    vA = mX \ vB;
    mA = mX \ mB;

    mA = mA .+ vA;

    return mA;
end

function LeastSquares( matrixSize, mX, mY  )

    mB = randn(matrixSize, matrixSize);
    vB = randn(matrixSize);

    mXT=transpose(mX);
    vA = ( mXT * mX) \ ( mXT * vB);
    mA = ( mXT * mX) \ ( mXT * mB);

    mA = mA .+ vA;

    return mA;
end

function CalcDistanceMatrix( matrixSize, mX, mY  )

    mY = randn(matrixSize, matrixSize);

    mA = transpose( sum(mX .^ 2, dims=1) ) .- (2 .* transpose(mX) * mY) .+ sum(mY .^ 2, dims=1);

    return mA;
end

function KMeans( matrixSize, mX, mY  )

    # Assuming Samples are slong Columns (Rows are features)
    numClusters     = Int64(max(round(matrixSize / 100), 1));
    vClusterId      = zeros(matrixSize);
    numIterations   = 100;

    # http://stackoverflow.com/questions/36047516/julia-generating-unique-random-integer-array
    mA          = mX[:, randperm(matrixSize)[1:numClusters]]; #<! Cluster Centroids

    for ii = 1:numIterations
        vMinDist, vClusterId[:] = findmin(transpose(sum(mA .^ 2, dims=1)) .- (2 .* transpose(mA)* mX), dims=1); #<! Is there a `~` equivalent in Julia?
        for jj = 1:numClusters
            mA[:, jj] = sum(mX[:, vClusterId .== jj], dims=2) ./ matrixSize;
        end
    end

    mA = mA[:, 1] .+ transpose(mA[:, end]);

    return mA;
end
