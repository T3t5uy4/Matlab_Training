% clc
% clear
% close all
warning('off');
% Used to add the current working directory and all its subdirectories to MATLAB's search path
addpath(genpath(pwd));

% Set initial training parameters
% initial population size
searchAgentsNum = 30;
numOfRecord = 40;
% population dimension
dim = 30;
% Maximum number of iterations, recommended maximum training is 300,000
maxFes = 300000;
% Fold is recommended to be set to 30
fold = 30;

%% Algorithm Select
% Select the algorithm to be trained
% Please do not select more than 13 algorithms
algorithmName = {'DEAHHO', 'DE', 'HHO', 'WOA', 'SCA', 'PSO', 'ABC', 'ALCPSO', 'SHADE', 'CLPSO', 'SADE', 'LSHADE'};
remarkStr = 'versionpk';

%% Select training data set
% 1-23 is the CEC05 function set
% 24-51 is the CEC13 function set
% 52-81 is the CEC14 function set
% 82-111 is the CEC17 function set and F83 has been delete
% 112-141 is the CEC 18 function set and F113 has been delete
% 142-151 is the CEC19 function set
% CEC2005
dataSetName = 'CEC05';
functionNameList = {'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'F17', 'F18', 'F19', 'F20', 'F21', 'F22', 'F23'};
% CEC2013
% dataSetName = 'CEC13';
% functionNameList = {'F24', 'F25', 'F26', 'F27', 'F28', 'F29', 'F30', 'F31', 'F32', 'F33', 'F34', 'F35', 'F36', 'F37', 'F38', 'F39', 'F40', 'F41', 'F42', 'F43', 'F44', 'F45', 'F46', 'F47', 'F48', 'F49', 'F50', 'F51'};
% CEC2014
% dataSetName = 'CEC14';
% functionNameList = {'F52', 'F53', 'F54', 'F55', 'F56', 'F57', 'F58', 'F59', 'F60', 'F61', 'F62', 'F63', 'F64', 'F65', 'F66', 'F67', 'F68', 'F69', 'F70', 'F71', 'F72', 'F73', 'F74', 'F75', 'F76', 'F77', 'F78', 'F79', 'F80', 'F81'};
% CEC2017
% dataSetName = 'CEC17';
% functionNameList = {'F82', 'F84', 'F85', 'F86', 'F87', 'F88', 'F89', 'F90', 'F91', 'F92', 'F93', 'F94', 'F95', 'F96', 'F97', 'F98', 'F99', 'F100', 'F101', 'F102', 'F103', 'F104', 'F105', 'F106', 'F107', 'F108', 'F109', 'F110', 'F111'};
% CEC 2018
% dataSetName = 'CEC18';
% functionNameList = {'F112', 'F114', 'F115', 'F116', 'F117', 'F118', 'F119', 'F120', 'F121', 'F122', 'F123', 'F124', 'F125', 'F126', 'F127', 'F128', 'F129', 'F130', 'F131', 'F132', 'F133', 'F134', 'F135', 'F136', 'F137', 'F138', 'F139', 'F140', 'F141'};
% CEC2019
% dataSetName = 'CEC19';
% functionNameList = {'F142', 'F143', 'F144', 'F145', 'F146', 'F147', 'F148', 'F149', 'F150', 'F151'};

isOutput = false;

% Prepare output preprocessing
% Get the current year, month, day, time and minute
dateStr = datestr(now, 'yyyy-mm-dd');
timeStr = datestr(now, 'HH_MM_SS');
% Set output directory name and file name
algNameStr = algorithmName{1};
dirName = ['output/', dateStr, '/', algNameStr, '_', remarkStr, '/', dataSetName, '_', timeStr, '_D=', num2str(dim), '_T=', num2str(maxFes), '_F=', num2str(fold)];
mkdir(dirName);
fileName = [dirName, '/', algNameStr];
% Excel header definition
te = {'F', 'algorithm', 'max', 'min', 'mean', 'std'};
xlsFileName = [fileName, '.xlsx'];
xlswrite(xlsFileName, te, 'overall')
% the startLineNum of overall sheet to write data
startLineNum = 2;

algorithmNum = size(algorithmName, 2);
% Generate the space of lineStyles, markerEdgeColors,Markers
nLines = algorithmNum;
basicLinestyles = cellstr(char('-', ':', '-.', '--'));
basicMarkers = cellstr(char('o', 'x', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p', 'h', '.'));
markerEdgeColors = hsv(nLines);
lineStyles = repmat(basicLinestyles, ceil(nLines / numel(basicLinestyles)), 1);
markers = repmat(basicMarkers, ceil(nLines / numel(basicMarkers)), 1);

for functionNum = 1:size(functionNameList, 2)
    functionName = functionNameList{functionNum};
    [lb, ub, dim, fobj] = getFunctions(functionName, dim);
    disp(['----------------', functionName, '----------------']);
    functionName = ['F', num2str(functionNum)];
    % Benchmark function
    resultCurves = zeros(algorithmNum, fold, numOfRecord);

    parfor cfold = 1:fold

        tic

        for cnum = 1:algorithmNum
            disp(['fold', num2str(cfold), '_', dataSetName, ':The ', algorithmName{cnum}, ' algorithm is trainning...']);
            algorithm = str2func(algorithmName{cnum});
            [~, ~, curve] = algorithm(searchAgentsNum, maxFes, lb, ub, dim, fobj);
            resultCurves(cnum, cfold, :) = uniformSampling(curve, numOfRecord);
        end

        toc

    end

    % Write data to excel sheet
    algorithmNamLabels = cell(algorithmNum * fold, 1);
    allCurves = zeros(algorithmNum * fold, numOfRecord);
    foldLables = repmat(int32(1:fold)', [algorithmNum, 1]);

    for it = 1:algorithmNum
        algorithmNamLabels((it - 1) * fold + 1:(it - 1) * fold + fold) = algorithmName(it);
        allCurves((it - 1) * fold + 1:(it - 1) * fold + fold, :) = resultCurves(it, :, :);
    end

    xlswrite(xlsFileName, algorithmNamLabels, functionName, 'A1')
    xlswrite(xlsFileName, foldLables, functionName, 'B1')
    xlswrite(xlsFileName, allCurves, functionName, 'C1')

    statisticValues = zeros(algorithmNum, 4);

    for it = 1:algorithmNum
        statisticValues(it, :) = [max(resultCurves(it, :, end)), min(resultCurves(it, :, end)), mean(resultCurves(it, :, end)), std(resultCurves(it, :, end))];
    end

    functionNumLable = repmat({functionName}, algorithmNum, 1);
    xlswrite(xlsFileName, functionNumLable, 'overall', ['A', num2str(startLineNum)])
    xlswrite(xlsFileName, algorithmName', 'overall', ['B', num2str(startLineNum)])
    xlswrite(xlsFileName, statisticValues, 'overall', ['C', num2str(startLineNum)])
    startLineNum = startLineNum + algorithmNum;

    % plot curveline
    clf
    set(gcf, 'Position', [0, 0, 1000, 600])

    for it = 1:algorithmNum
        yy(it, :) = mean(allCurves((it - 1) * fold + 1:(it - 1) * fold + fold, :));
    end

    xx = [1:numOfRecord] * (maxFes / numOfRecord);

    for it = 1:algorithmNum
        semilogy(xx, yy(it, :), [lineStyles{it} markers{it}], 'LineWidth', 1.5, 'Color', markerEdgeColors(it, :));
        hold on;
    end

    hold off;
    title(functionName);
    set(gcf, 'color', 'white')
    set(gca, 'YScale', 'log', 'YLimMode', 'auto')
    xlabel('t');
    ylabel('Best fitness');
    algorithmName1 = strrep(algorithmName, '_', '\_');

    legend(algorithmName1, 'Location', 'northeast');
    legend('boxoff')

    a = findobj(gcf); % get the handles associated with the current figure
    allaxes = findall(a, 'Type', 'axes');
    % alllines=findall(a,'Type','line');
    alltext = findall(a, 'Type', 'text');
    set(allaxes, 'FontName', 'Times', 'LineWidth', 1, 'FontSize', 14, 'FontWeight', 'bold');
    % set(alllines,'Linewidth',1);
    set(alltext, 'FontName', 'Times', 'FontSize', 14, 'FontWeight', 'bold')
    %
    set(gcf, 'PaperUnits', 'inches');
    % set(gcf, 'PaperUnits', 'centimeters');
    krare = 3.5;
    % x_width=krare*1.618 ;
    x_width = krare * 5/3;
    %  x_width=3*1;
    y_width = krare * 4/3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 4])
    set(gca, ...
        'Box', 'on', ...
        'TickDir', 'in', ...
        'TickLength', [.02 .02], ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'YGrid', 'off', ...
        'XGrid', 'off', ...
        'XColor', [.3 .3 .3], ...
        'YColor', [.3 .3 .3], ...
        'LineWidth', 1);
    axis tight
    %     grid on
    %     box on
    saveas(gcf, [fileName, '-', functionName, '-carve'], 'fig')
    fileName1 = [fileName, '-', functionName, '-carve'];
    print(fileName1, '-dtiff', '-r300'); %<-Save as PNG with 300 DPI

end

Orderhao(xlsFileName);
pValueToExcelhao(xlsFileName, fold);
FridTest3(xlsFileName, fold)
% FridTest4(xlsFileName, fold)
disp([dataSetName, ':', algNameStr, ' training is over!'])

close all
