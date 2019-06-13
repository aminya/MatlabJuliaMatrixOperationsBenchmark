function [ mRunTime ] = MatlabBench( operationMode, vTestIdx )

FALSE   = 0;
TRUE    = 1;
OFF     = 0;
ON      = 1;

RUN_TIME_DATA_FOLDER    = 'RunTimeData';
RUN_TIME_FILE_NAME      = 'RunTimeMatlab.csv';

OPERATION_MODE_PARTIAL  = 1; %<! For Testing (Runs Fast)
OPERATION_MODE_FULL     = 2;

cRunTimeFunctionsBase = {@MatrixGeneration, @MatrixAddition,@MatrixMultiplication,...
    @MatrixQuadraticForm, @MatrixReductions, @ElementWiseOperations, @MatrixExp,...
    @MatrixSqrt, @Svd, @Eig,@CholDec, @MatInv, @LinearSystem, @LeastSquares,...
    @CalcDistanceMatrix, @KMeans};

cFunctionStringBase = {'Matrix Generation', 'Matrix Addition', 'Matrix Multiplication',...
    'Matrix Quadratic Form', 'Matrix Reductions', 'Element Wise Operations',...
    'Matrix Exponential', 'Matrix Square Root', 'SVD', 'Eigen Decomposition',...
    'Cholesky Decomposition', 'Matrix Inversion','Linear System Solution',...
    'Linear Least Squares', 'Squared Distance Matrix', 'K-Means Run Time'};

numTests = length(cRunTimeFunctionsBase);

if(operationMode == OPERATION_MODE_PARTIAL)
    vMatrixSize = csvread('vMatrixSizePartial.csv');
    numIterations = csvread('numIterationsPartial.csv');
elseif(operationMode == OPERATION_MODE_FULL)
    vMatrixSize = csvread('vMatrixSizeFull.csv');
    numIterations = csvread('numIterationsFull.csv');
end

cRunTimeFunctions = cRunTimeFunctionsBase(vTestIdx);
cFunctionString = cFunctionStringBase(vTestIdx);

mRunTime = zeros(length(vMatrixSize), length(cRunTimeFunctions), numIterations);

for ii = 1:length(vMatrixSize)
    matrixSize = vMatrixSize(ii);
    mX = randn(matrixSize, matrixSize);
    mY = randn(matrixSize, matrixSize);
    disp(['Matrix Size - ', num2str(matrixSize)]);
    
    for jj = 1:length(cRunTimeFunctions)
        disp(['Processing ', num2str(cFunctionString{jj}), ' Matrix Size ', num2str(matrixSize)]);
        
        for kk = 1:numIterations
            
            fun=@() cRunTimeFunctions{jj}(matrixSize, mX, mY);
            mRunTime(ii, jj, kk) = timeit(fun); % computes median of bench times
        
        end
    end
end

runTimeFilePath = fullfile(RUN_TIME_DATA_FOLDER, RUN_TIME_FILE_NAME);

mRunTimeBase = 0;
if(exist(runTimeFilePath, 'file'))
    mRunTimeBase = csvread(runTimeFilePath);
end

if(any(size(mRunTimeBase) ~= [length(vMatrixSize), numTests]))
    % Previous Data has incompatible dimensions
    mRunTimeBase = zeros([length(vMatrixSize), numTests]);
end

mRunTimeBase(:, vTestIdx) = mRunTime;

if(operationMode == OPERATION_MODE_FULL)
    csvwrite(runTimeFilePath, mRunTimeBase);
end


end


function mA = MatrixGeneration( matrixSize, mX, mY )

mA = randn(matrixSize, matrixSize);
mB = rand(matrixSize, matrixSize);

end

function mA = MatrixAddition( matrixSize, mX, mY )

scalarA = rand(1);
scalarB = rand(1);

mA = (scalarA .* mX) + (scalarB .* mY);

end

function mA = MatrixMultiplication( matrixSize, mX, mY )

sacalrA = rand(1);
sacalrB = rand(1);

mA = (sacalrA + mX) * (sacalrB + mY);

end

function mA = MatrixQuadraticForm( matrixSize, mX, mY )

vX = randn(matrixSize, 1);
vB = randn(matrixSize, 1);
sacalrC = rand(1);

mA = ((mX * vX).' * (mX * vX)) + (vB.' * vX) + sacalrC;

end

function mA = MatrixReductions( matrixSize, mX, mY )

mA = sum(mX, 1) + min(mY, [], 2);

end

function mA = ElementWiseOperations( matrixSize, mX, mY )

% Make sure roots are positive
mA = rand(matrixSize, matrixSize);
mB = 3 + rand(matrixSize, matrixSize);
mC = rand(matrixSize, matrixSize);

mD = abs(mA) + sin(mB);
mE = exp(-(mA .^ 2));
mF = (-mB + sqrt((mB .^ 2) - (4 .* mA .* mC))) ./ (2 .* mA);

mA = mD + mE + mF;

end
function mA = MatrixExp(matrixSize, mX, mY)

mA = expm(mX);

end

function mA = MatrixSqrt(matrixSize, mX, mY)

mY = mX.' * mX;

mA = sqrtm(mY);

end

function mA = Svd(matrixSize, mX, mY) % check

[mU, mS, mV] = svd(mX);

mA = mU + mS + mV;

end

function mA = Eig(matrixSize, mX, mY) % check

[mD, mV] = eig(mX);

mA = mD + mV;

end

function mA = CholDec(matrixSize, mX, mY)

mY = mX.' * mX;

mA = chol(mY);
  
end

function mA = MatInv(matrixSize, mX, mY)

mY = mX.' * mX;

mA = inv(mY);
mB = pinv(mX);

mA = mA + mB;

end


function mA = LinearSystem(matrixSize, mX, mY)

% mA = randn(matrixSize, matrixSize);
mB = randn(matrixSize, matrixSize);
vB = randn(matrixSize, 1);

vA = mX \ vB;
mA = mX \ mB;

mA = mA + vA;

end

function mA = LeastSquares(matrixSize, mX, mY)

% mA = randn(matrixSize, matrixSize);
mB = randn(matrixSize, matrixSize);
vB = randn(matrixSize, 1);

vA = (mX.' * mX) \ (mX.' * vB);
mA = (mX.' * mX) \ (mX.' * mB);

mA = mA + vA;

end

function mA = CalcDistanceMatrix(matrixSize, mX, mY)

mY = randn(matrixSize, matrixSize);

mA = sum(mX .^ 2, 1).' - (2 .* mX.' * mY) + sum(mY .^ 2, 1);

end

function mA = KMeans(matrixSize, mX, mY)

% Assuming Samples are slong Columns (Rows are features)
% mX              = randn(matrixSize, matrixSize);
numClusters     = max(round(matrixSize / 100), 1);
vClusterId      = zeros(matrixSize, 1);
numIterations   = 10;


mA  = mX(:, randperm(matrixSize, numClusters)); %<! Cluster Centroids

for ii = 1:numIterations
    [~, vClusterId(:)] = min(sum(mA .^ 2, 1).' - (2 .* mA.' * mX), [], 1);
    for jj = 1:numClusters
        mA(:, jj) = sum(mX(:, vClusterId == jj), 2) ./ sum(vClusterId == jj);
    end
    
end

mA = mA(:, 1) + mA(:, end).';
  
end


