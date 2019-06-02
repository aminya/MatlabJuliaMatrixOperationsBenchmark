# ----------------------------------------------------------------------------------------------- #
# Julia Matrix Operations Benchmark - Test Suite 0003
# Reference:
#   1. C.
# Remarks:
#   1.  W.
# TODO:
#   1.  A
#   Release Notes:
#   -   1.0.002     10/02/2017  Royi Avital
#       *   Added generation of 'mX' once outside the functions.
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

function JuliaMatrixBenchmark0003( operationMode = 2 )

  OPERATION_MODE_PARTIAL  = 1; # For Testing (Runs Fast)
  OPERATION_MODE_FULL     = 2;

  cRunTimeFunctions = [LinearSystemRunTime, LeastSquaresRunTime, CalcDistanceMatrixRunTime, KMeansRunTime];

  cFunctionString = ["Linear System Solution", "Linear Least Squares", "Squared Distance Matrix", "K-Means"]

  if(operationMode == OPERATION_MODE_PARTIAL)
    vMatrixSize = round(Int64, squeeze(readcsv("vMatrixSizePartial.csv"), 1));
    numIterations = round(Int64, squeeze(readcsv("numIterationsPartial.csv"), 1));
  elseif(operationMode == OPERATION_MODE_FULL)
    vMatrixSize = round(Int64, squeeze(readcsv("vMatrixSizeFull.csv"), 1));
    numIterations = round(Int64, squeeze(readcsv("numIterationsFull.csv"), 1));
  end

  numIterations = numIterations[1]; # It is 1x1 Array -> Scalar

  mRunTime = zeros(length(vMatrixSize), length(cRunTimeFunctions), numIterations);

  tic();
  for ii = 1:length(vMatrixSize)
    matrixSize = vMatrixSize[ii];
    mX = randn(matrixSize, matrixSize);
    println("Matrix Size - $matrixSize");
    for jj = 1:length(cRunTimeFunctions)
      println("Processing $(cFunctionString[jj]) Matrix Size $matrixSize");
      for kk = 1:numIterations
        mA, mRunTime[ii, jj, kk] = cRunTimeFunctions[jj](matrixSize, mX);
      end
      println("Finished Processing $(cFunctionString[jj])");
    end
  end
  totalRunTime = toq();

  mRunTime = median(mRunTime, 3);
  mRunTime = squeeze(mRunTime, 3);

  println("Finished the Benchmark in $totalRunTime [Sec]");

  writecsv("RunTimeJulia0003.csv", mRunTime);

  return mRunTime;

end

function LinearSystemRunTime( matrixSize, mX )

  mB = randn(matrixSize, matrixSize);
  vB = randn(matrixSize);

  tic();
  vA = mX \ vB;
  mA = mX \ mB;
  runTime = toq();

  mA = mA .+ vA;

  return mA, runTime;
end

function LeastSquaresRunTime( matrixSize, mX )

  mB = randn(matrixSize, matrixSize);
  vB = randn(matrixSize);

  tic();
  mXT=transpose(mX);
  vA = ( mXT * mX) \ ( mXT * vB);
  mA = ( mXT * mX) \ ( mXT * mB);
  runTime = toq();

  mA = mA .+ vA;

  return mA, runTime;
end

function CalcDistanceMatrixRunTime( matrixSize, mX )

  mY = randn(matrixSize, matrixSize);

  tic();
  mA = transpose(sum(mX .^ 2, 1)) .- (2 .* transpose(mX) * mY) .+ sum(mY .^ 2, 1);
  runTime = toq();

  return mA, runTime;
end

function KMeansRunTime( matrixSize, mX )

  # Assuming Samples are slong Columns (Rows are features)
  numClusters     = Int64(max(round(matrixSize / 100), 1));
  vClusterId      = zeros(matrixSize);
  numIterations   = 100;

  tic();
  # http://stackoverflow.com/questions/36047516/julia-generating-unique-random-integer-array
  mA          = mX[:, randperm(matrixSize)[1:numClusters]]; #<! Cluster Centroids

  for ii = 1:numIterations
    vMinDist, vClusterId[:] = findmin(transpose(sum(mA .^ 2, 1)) .- (2 .* transpose(mA)* mX), 1); #<! Is there a `~` equivalent in Julia?
    for jj = 1:numClusters
      mA[:, jj] = sum(mX[:, vClusterId .== jj], 2) ./ matrixSize;
    end
  end

  runTime = toq();

  mA = mA[:, 1] .+ transpose(mA[:, end]);

  return mA;
end
