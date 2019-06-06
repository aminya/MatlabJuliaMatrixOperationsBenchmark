

function JuliaMatrixBenchmarkOpt0001( operationMode = 2 )

  allFunctions = [MatrixGeneration, MatrixAddition, MatrixMultiplication, MatrixQuadraticForm, MatrixReductions, ElementWiseOperations, MatrixExp, MatrixSqrt, Svd, Eig, CholDec, MatInv, LinearSystem, LeastSquares];

  # need tweaking:
  # , CalcDistanceMatrix, KMeans

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

        benchijk =@benchmark $fun($matrixSize, $mX, $mY);
        mRunTime[ii, jj, kk]=minimum(benchijk.times);

      end
      jj+=1;
    end

  end

  writedlm("RunTimeJuliaOptimized001.csv", mRunTime);

  return mRunTime;
end

function MatrixGeneration( matrixSize, mX, mY )

  mA = randn(matrixSize, matrixSize);
  mB = rand(matrixSize, matrixSize);

  mA = mA .+ mB;

  return mA;
end

function MatrixAddition( matrixSize, mX, mY )

  scalarA = rand();
  scalarB = rand();

  # mA = (scalarA .* mX) .+ (scalarB .* mY);
  mA = Array(Float64, matrixSize, matrixSize);
  @simd for ii = 1:(matrixSize * matrixSize)
    @inbounds mA[ii] = (scalarA * mX[ii]) + (scalarB * mY[ii]);
  end

  return mA;
end

function MatrixMultiplication( matrixSize, mX, mY )

  scalarA = rand();
  scalarB = rand();

  # mA = (scalarA .+ mX) * (scalarB .+ mY);
  mA = Array(Float64, matrixSize, matrixSize);
  @simd for ii = 1:(matrixSize * matrixSize)
    @inbounds mA[ii] = (scalarA + mX[ii]) * (scalarB + mY[ii]);
  end

  return mA;
end

function MatrixQuadraticForm( matrixSize, mX, mY )

  vX = randn(matrixSize);
  vB = randn(matrixSize);
  sacalrC = rand();

  mA = (transpose(mX * vX) * (mX * vX)) .+ (transpose(vB) * vX) .+ sacalrC;

  return mA;
end

function MatrixReductions( matrixSize, mX, mY )

  mA = sum(mX, 1) .+ minimum(mY, 2); #Broadcasting

  return mA;
end

function ElementWiseOperations( matrixSize, mX, mY )

  mA = rand(matrixSize, matrixSize);
  mB = 3 .+ rand(matrixSize, matrixSize);
  mC = rand(matrixSize, matrixSize);

  # mD = abs.(mA) .+ sin.(mA);
  mD = Array(Float64, matrixSize, matrixSize);
  @simd for ii = 1:(matrixSize * matrixSize)
    @inbounds mD[ii] = abs(mA[ii]) + sin(mA[ii]);
  end

  # mE = exp.(-(mA .^ 2));
  mE = Array(Float64, matrixSize, matrixSize);
  @simd for ii = 1:(matrixSize * matrixSize)
    @inbounds mE[ii] = exp(- (mA[ii] * mA[ii]));
  end

  # mF = (-mB .+ sqrt.((mB .^ 2) .- (4 .* mA .* mC))) ./ (2 .* mA);
  mF = Array(Float64, matrixSize, matrixSize);
  @simd for ii = 1:(matrixSize * matrixSize)
    @inbounds mF[ii] = (-mB[ii] + sqrt( (mB[ii] * mB[ii]) - (4 * mA[ii] * mC[ii]) )) ./ (2 * mA[ii]);
  end

  mA = mD .+ mE .+ mF;

  return mA;
end


function MatrixExp( matrixSize, mX, mY )

  mA = expm(mX);

  return mA;
end

function MatrixSqrt( matrixSize, mX, mY )

  mY = transpose(mX) * mX;

  mA = sqrtm(mY);

  return mA;
end

function Svd( matrixSize, mX, mY )

  mU, mS, mV = svd(mX, thin = false);

  mA = mU .+ mS .+ mV;

  return mA;
end

function Eig( matrixSize, mX, mY )

  mD, mV = eig(mX);

  mA = mD .+ mV;

  return mA;
end

function CholDec( matrixSize, mX, mY )

  mY = transpose(mX) * mX;

  mA = chol(mY);

  return mA;
end

function MatInv( matrixSize, mX, mY )

  mY = transpose(mX) * mX;

  mA = inv(mY);
  mB = pinv(mX);

  mA = mA .+ mB;

  return mA;
end

function LinearSystem( matrixSize, mX, mY )

  mB = randn(matrixSize, matrixSize);
  vB = randn(matrixSize);

  vA = mX \ vB;
  mA = mX \ mB;

  mA = mA .+ vA;

  return mA;
end

function LeastSquares( matrixSize, mX, mY )

  mB = randn(matrixSize, matrixSize);
  vB = randn(matrixSize);

  vA = (transpose(mX) * mX) \ (transpose(mX) * vB);
  mA = (transpose(mX) * mX) \ (transpose(mX) * mB);

  mA = mA .+ vA;

  return mA;
end

function CalcDistanceMatrix( matrixSize, mX, mY )

  mY = randn(matrixSize, matrixSize);

  mA = transpose(sum(mX .^ 2, 1)) .- (2 .* transpose(mX) * mY) .+ sum(mY .^ 2, 1);

  return mA;
end

function KMeans( matrixSize, mX, mY )

  # Assuming Samples are slong Columns (Rows are features)
  numClusters     = Int64(max(round(matrixSize / 100), 1));
  vClusterId      = zeros(matrixSize);
  numIterations   = 10;

  # http://stackoverflow.com/questions/36047516/julia-generating-unique-random-integer-array
  mA          = mX[:, randperm(matrixSize)[1:numClusters]]; #<! Cluster Centroids

  for ii = 1:numIterations
    vMinDist, vClusterId[:] = findmin(sum(mA .^ 2, 1).' .- (2 .* transpose(mA) * mX), 1); #<! Is there a `~` equivalent in Julia?
    for jj = 1:numClusters
      mA[:, jj] = sum(mX[:, vClusterId .== jj], 2) ./ sum(vClusterId .== jj);
    end
  end


  mA = mA[:, 1] .+ transpose(mA[:, end]);

  return mA;
end
