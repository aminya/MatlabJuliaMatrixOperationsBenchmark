% ----------------------------------------------------------------------------------------------- %
% Analyze MATLAB & Julia Run Time Results
% Reference:
%   1. C.
% Remarks:
%   1.  W.
% TODO:
%   1.  A
%   Release Notes:
%   -   1.0.001     11/02/2017  Royi Avital
%       *   Added support for Julia Optimized.
%       *   Saving figures into Figures sub folder.
%   -   1.0.000     09/02/2017  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

%% Setting Enviorment Parameters

run('InitScript.m');

RUN_TIME_DATA_FOLDER = 'RunTimeData\';

MATLAB_IDX  = 1;
JULIA_IDX   = 2;

MATLAB_RUN_TIME_FILE_NAME       = 'RunTimeMatlab.csv';
JULIA_RUN_TIME_FILE_NAME        = 'RunTimeJulia.csv';
JULIA_OPT_RUN_TIME_FILE_NAME    = 'RunTimeJuliaOpt.csv';

cLegendString = {'MATLAB', 'Julia', 'Julia Optimized'};

figureIdx           = 0;
figureCounterSpec   = '%04d';

vMatrixSize = csvread('vMatrixSizeFull.csv');
cFunctionString = {'Matrix Generation','Matrix Addition','Matrix Multiplication',...
    'Matrix Quadratic Form', 'Matrix Reductions', 'Element Wise Operations',...
    'Matrix Exponential', 'Matrix Square Root', 'SVD','Eigen Decomposition',...
    'Cholesky Decomposition', 'Matrix Inversion','Linear System Solution',...
    'Linear Least Squares', 'Squared Distance Matrix', 'K-Means'};



%% Setting Parameters

generateImages = ON;


%% Loading Data

mRunTimeMatlab      = csvread([RUN_TIME_DATA_FOLDER, MATLAB_RUN_TIME_FILE_NAME]);
mRunTimeJulia       = csvread([RUN_TIME_DATA_FOLDER, JULIA_RUN_TIME_FILE_NAME]);
mRunTimeJuliaOpt    = csvread([RUN_TIME_DATA_FOLDER, JULIA_OPT_RUN_TIME_FILE_NAME]);

numTests    = length(cFunctionString);
numMatSize  = length(vMatrixSize);

if(any(size(mRunTimeMatlab) ~= size(mRunTimeJulia)))
    error('Run Time Data Dimensions Don''t Match');
end

if(any(size(mRunTimeMatlab) ~= size(mRunTimeJuliaOpt)))
    error('Run Time Data Dimensions Don''t Match');
end

if(size(mRunTimeMatlab, 2) ~= numTests)
    error('Run Time Data Has Incompatible Number of Tests');
end

if(size(mRunTimeMatlab, 1) ~= numMatSize)
    error('Run Time Data Has Incompatible Number of Matrix Size');
end


%% Displaying Results

for ii = 1:numTests
    
    figureIdx   = figureIdx + 1;
    hFigure     = figure('Position', figPosMedium);
    hAxes       = axes();
    set(hAxes, 'NextPlot', 'add');
    hLineSeries = plot(vMatrixSize, [mRunTimeMatlab(:, ii), mRunTimeJulia(:, ii), mRunTimeJuliaOpt(:, ii)]);
    set(hLineSeries, 'LineWidth', lineWidthNormal);
    set(get(hAxes, 'Title'), 'String', ['Test - ', cFunctionString{ii}], ...
        'FontSize', fontSizeTitle);
    set(get(hAxes, 'XLabel'), 'String', 'Matrix Size', ...
        'FontSize', fontSizeAxis);
    set(get(hAxes, 'YLabel'), 'String', 'Run Time  [Sec]', ...
        'FontSize', fontSizeAxis);
    hLegend = ClickableLegend(cLegendString);
    
    if(generateImages == ON)
        set(hAxes, 'LooseInset', [0.05, 0.05, 0.05, 0.05]);
        saveas(hFigure,['Figures\Figure', num2str(figureIdx, figureCounterSpec), '.png']);
    end

end


%% Restoring Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLooseInset);

