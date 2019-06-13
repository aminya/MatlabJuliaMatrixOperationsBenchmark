

function JuliaMatrixBenchmarkOpt0001( operationMode = 2 )

  allFunctions = [MatrixGeneration, MatrixAddition, MatrixMultiplication, MatrixQuadraticForm, MatrixReductions, ElementWiseOperations, MatrixExp, MatrixSqrt, Svd, Eig, CholDec, MatInv, LinearSystem, LeastSquares, CalcDistanceMatrix, KMeans];


  allFunctionsString = ["Matrix Generation", "Matrix Addition", "Matrix Multiplication", "Matrix Quadratic Form", "Matrix Reductions", "Element Wise Operations", "Matrix Exponential", "Matrix Square Root", "Singular Value Decomposition", "Eigen Decomposition","Cholesky Decomposition", "Matrix Inversion", "Linear System Solution", "Linear Least Squares", "Squared Distance Matrix", "K-Means"];

  if (operationMode == 1) # partial benchmark
    vMatrixSize =  dropdims(DelimitedFiles.readdlm("vMatrixSizePartial.csv", ',',Int64), dims=1);
    numIterations = dropdims(DelimitedFiles.readdlm("numIterationsPartial.csv", ',',Int64), dims=1);

  elseif (operationMode == 2) # full benchmark
    vMatrixSize = dropdims(readdlm("vMatrixSizeFull.csv", ',',Int64), dims=1);
    numIterations =  dropdims(readdlm("numIterationsFull.csv", ',',Int64), dims=1);

  elseif (operationMode == 0) # Test benchmark
    vMatrixSize = 2;
    numIterations =  1;

  end

  numIterations = numIterations[1]; # It is 1x1 Array -> Scalar

  mRunTime = zeros(length(vMatrixSize), length(allFunctions), numIterations);
  tRunTime= Array{Any}(undef,length(mRunTime)+1,3)# a table containing all the information
  tRunTime[1,:]=["Function Name","Matrix Size","Run Time"];

  rr=2; # row counter for table

  for ii = 1:length(vMatrixSize)
    matrixSize = vMatrixSize[ii];
    mX = randn(matrixSize, matrixSize);
    mY = randn(matrixSize, matrixSize);
    println("Matrix Size - $matrixSize");

    jj=1;
    for fun in allFunctions
      println("Processing $(allFunctionsString[jj]) - MatrixSize= $matrixSize");

      for kk = 1:numIterations;

        benchIJK =@benchmark $fun($matrixSize, $mX, $mY);

        mRunTime[ii, jj, kk]=median(benchIJK.times);
        tRunTime[rr,:]=["$(allFunctionsString[jj])","$matrixSize",mRunTime[ii, jj, kk]];

        rr+=1;
      end
      jj+=1;
    end

  end
  writedlm("RunTimeJuliaTable.csv", tRunTime,',');
  writedlm("RunTimeJulia.csv", mRunTime,',');

  return tRunTime;
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

  mA = sum(mX, dims=1) .+ minimum(mY, dims=2); #Broadcasting

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

  F = svd(mX, full = false); # F is SVD object
  mU, mS, mV = F;

  return mA=0;
end

function Eig( matrixSize, mX, mY )

  F  = eigen(mX); # F is eigen object
  mD, mV = F;

  return mA=0;
end

function CholDec( matrixSize, mX, mY )

  mY = transpose(mX) * mX;

  mA = cholesky(mY);

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

  mA = transpose(sum(mX .^ 2, dims=1)) .- (2 .* transpose(mX) * mY) .+ sum(mY .^ 2, dims=1);

  return mA;
end

function KMeans( matrixSize, mX, mY )

  # Assuming Samples are slong Columns (Rows are features)
  numClusters     = Int64( max( round(matrixSize / 100), 1 ) ); # % max between 1 and round(...)
  numIterations   = 10;

  # http://stackoverflow.com/questions/36047516/julia-generating-unique-random-integer-array
  mA          = mX[:, randperm(matrixSize)[1:numClusters]]; #<! Cluster Centroids

  for ii = 1:numIterations
    vMinDist, vClusterId = findmin( transpose(sum(mA .^ 2, 1)) .- (2 .* transpose(mA) * mX), dims=1); #<! Is there a `~` equivalent in Julia?
    vClusterId=dropdims(mClusterId, dims=1); # to be able to access it later
    for jj = 1:numClusters
      mA[:, jj] = sum(mX[:, vClusterId .== jj], dims=2) ./ sum(vClusterId .== jj);
    end
  end


  mA = mA[:, 1] .+ transpose(mA[:, end]);

  return mA;
end
