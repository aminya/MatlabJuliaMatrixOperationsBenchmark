function [ mRunTime ] = MatlabMatrixBenchmark0001( operationMode, vTestIdx )
% ----------------------------------------------------------------------------------------------- %
% MATLAB Matrix Operations Benchmark - Test Suite 0001
% Reference:
%   1. C.
% Remarks:
%   1.  Keep 'mX' and 'mY' "Read Only" within the functions to match Julia (Pass by Address).
% TODO:
%   1.  A
%   Release Notes:
%	- 	1.0.003 	12/02/2017	Royi Avital
% 		* 	Ability to run only some of the tests.
%   -   1.0.002     10/02/2017  Royi Avital
%       *   Added generation of 'mX' and 'mY' once outside the functions.
%       *   Fixed issue with the Quadratic Form.
%   -   1.0.001     09/02/2017  Royi Avital
%       *   Added 'MatrixExpRunTime()' and 'MatrixSqrtRunTime()'.
%       *   Added Quadratic Matrix Form Calculation 'MatrixQuadraticFormRunTime()'.
%       *   Added Univariate Quadratic Function Root to 'ElementWiseOperationsRunTime()'.
%       *   Updated 'MatrixGenerationRunTime()' to include Uniform Random Number Generation.
%       *   Fixed issue with 'CalcDistanceMatrixRunTime'.
%   -   1.0.000     09/02/2017  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

FALSE   = 0;
TRUE    = 1;
OFF     = 0;
ON      = 1;

RUN_TIME_DATA_FOLDER    = 'RunTimeData';
RUN_TIME_FILE_NAME      = 'RunTimeMatlab0001.csv';

OPERATION_MODE_PARTIAL  = 1; %<! For Testing (Runs Fast)
OPERATION_MODE_FULL     = 2;

cRunTimeFunctionsBase = {@MatrixGenerationRunTime, @MatrixAdditionRunTime, @MatrixMultiplicationRunTime, ...
    @MatrixQuadraticFormRunTime, @MatrixReductionsRunTime, @ElementWiseOperationsRunTime};

cFunctionStringBase = {['Matrix Generation'], ['Matrix Addition'], ['Matrix Multiplication'], ...
    ['Matrix Quadratic Form'], ['Matrix Reductions'], ['Element Wise Operations']};

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

hTotelRunTimer = tic();
for ii = 1:length(vMatrixSize)
    matrixSize = vMatrixSize(ii);
    mX = randn(matrixSize, matrixSize);
    mY = randn(matrixSize, matrixSize);
    disp(['Matrix Size - ', num2str(matrixSize)]);
    for jj = 1:length(cRunTimeFunctions)
        disp(['Processing ', num2str(cFunctionString{jj}), ' Matrix Size ', num2str(matrixSize)]);
        for kk = 1:numIterations
            [mA, mRunTime(ii, jj, kk)] = cRunTimeFunctions{jj}(matrixSize, mX, mY);
        end
        disp(['Finished Processing ', num2str(cFunctionString{jj})]);
    end
end
totalRunTime = toc(hTotelRunTimer);

mRunTime = median(mRunTime, 3);

disp(['Finished the Benchmark in ', num2str(totalRunTime), ' [Sec]']);

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


function [ mA, runTime ] = MatrixGenerationRunTime( matrixSize, mX, mY )

tic();
mA = randn(matrixSize, matrixSize);
mB = rand(matrixSize, matrixSize);
runTime = toc();


end

function [ mA, runTime ] = MatrixAdditionRunTime( matrixSize, mX, mY )

scalarA = rand(1);
scalarB = rand(1);

tic();
mA = (scalarA .* mX) + (scalarB .* mY);
runTime = toc();


end

function [ mA, runTime ] = MatrixMultiplicationRunTime( matrixSize, mX, mY )

sacalrA = rand(1);
sacalrB = rand(1);

tic();
mA = (sacalrA + mX) * (sacalrB + mY);
runTime = toc();


end

function [ mA, runTime ] = MatrixQuadraticFormRunTime( matrixSize, mX, mY )

vX = randn(matrixSize, 1);
vB = randn(matrixSize, 1);
sacalrC = rand(1);

tic();
mA = ((mX * vX).' * (mX * vX)) + (vB.' * vX) + sacalrC;
runTime = toc();


end

function [ mA, runTime ] = MatrixReductionsRunTime( matrixSize, mX, mY )

tic();
mA = sum(mX, 1) + min(mY, [], 2);
runTime = toc();


end

function [ mA, runTime ] = ElementWiseOperationsRunTime( matrixSize, mX, mY )

% Make sure roots are positive
mA = rand(matrixSize, matrixSize);
mB = 3 + rand(matrixSize, matrixSize);
mC = rand(matrixSize, matrixSize);

tic();
mD = abs(mA) + sin(mB);
mE = exp(-(mA .^ 2));
mF = (-mB + sqrt((mB .^ 2) - (4 .* mA .* mC))) ./ (2 .* mA);
runTime = toc();

mA = mD + mE + mF;


end

