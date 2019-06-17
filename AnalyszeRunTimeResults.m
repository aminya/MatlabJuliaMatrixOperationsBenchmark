%% Setting Enviorment Parameters

AnalysisInitScript;

figureCounterSpec   = '%04d';

generateImages=1;

%% Loading Data

tRunTimeMatlab = readtable(fullfile('RunTimeData\', 'RunTimeMatlabTable.csv'));
mRunTimeMatlab=table2array(tRunTimeMatlab(2:end,2:end));
vMatrixSizeMatlab=table2array(tRunTimeMatlab(1,2:end));
sFunNameMatlab=table2array(tRunTimeMatlab(2:end,1));

tRunTimeJuliamkl = readtable(fullfile('RunTimeData\', 'RunTimeJuliamklTable.csv'));
mRunTimeJuliamkl=table2array(tRunTimeJuliamkl(2:end,2:end));
vMatrixSizeJuliamkl=table2array(tRunTimeJuliamkl(1,2:end));
sFunNameJuliamkl=table2array(tRunTimeJuliamkl(2:end,1));

tRunTimeJuliamklSIMD = readtable(fullfile('RunTimeData\', 'RunTimeJuliamklSIMDTable.csv'));
mRunTimeJuliamklSIMD=table2array(tRunTimeJuliamklSIMD(2:end,2:end));
vMatrixSizeJuliamklSIMD=table2array(tRunTimeJuliamklSIMD(1,2:end));
sFunNameJuliamklSIMD=table2array(tRunTimeJuliamklSIMD(2:end,1));

%% Displaying Results
figureIdx           = 0;

for ii = 1:size(mRunTimeMatlab,1)

    figureIdx   = figureIdx + 1;
    hFigure     = figure('Position', figPosMedium);
    hAxes       = axes();

    loglog(vMatrixSizeMatlab,mRunTimeMatlab(ii,:),'-o','LineWidth',lineWidthThin,'MarkerFaceColor','b'); hold on;
    loglog(vMatrixSizeJuliamkl,mRunTimeJuliamkl(ii,:),'-s','LineWidth',lineWidthThin,'MarkerFaceColor','r'); hold on;
    plotJuliaSIMD=ismember( sFunNameMatlab{ii}, sFunNameJuliamklSIMD ); % if 1 will plot JuliaSIMD
    if any(plotJuliaSIMD)
        h=loglog(vMatrixSizeJuliamklSIMD,mRunTimeJuliamklSIMD(plotJuliaSIMD,:),'-*','LineWidth',lineWidthThin,'MarkerFaceColor','y');
        legend('MATLAB','Julia-MKL','Julia-MKL-SIMD','Location','southeast')
    else
        legend('MATLAB','Julia-MKL','Location','southeast')
    end
    hold off;
    title(num2str(sFunNameMatlab{ii}));
    xlabel('Matrix Size');
    ylabel('Run Time  [micro Seconds]');

    if(generateImages == 1)
        set(hAxes, 'LooseInset', [0.05, 0.05, 0.05, 0.05]);
        saveas(hFigure,['Figures\Figure', num2str(figureIdx, figureCounterSpec), '.png']);
    end

end


%% Restoring Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLooseInset);
