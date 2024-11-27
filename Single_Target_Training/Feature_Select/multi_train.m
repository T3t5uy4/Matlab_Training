function multi_train(dataset)
    warning off;
    clc;
    close all;
    addpath(genpath(pwd));
    algorithms_str = {
                      'bDEAHHO', ...
                          'bDE', ...
                          'bHHO', ...
                          'bWOA', ...
                          'bGA', ...
                          'bBA', ...
                          'bPSO', ...
                          'bABC', ...
                          'bBMWOA', ...
                          'bCDLOBA', ...
                          'bOBSCA', ...
                          'bHHODE', ...
                          'bRLHHO', ...
                          'bALCPSO', ...
                          'bCLPSO', ...
                          'bSADE' ...
                      };
    bDEAHHO = 1;
    bDE = 1;
    bHHO = 1;
    bWOA = 0;
    bGA = 0;
    bBA = 0;
    bPSO = 0;
    bABC = 0;
    bBMWOA = 0;
    bCDLOBA = 0;
    bOBSCA = 0;
    bHHODE = 0;
    bRLHHO = 0;
    bALCPSO = 0;
    bCLPSO = 0;
    bSADE = 0;
    algorithms_colle = [
                        bDEAHHO, ...
                            bDE, ...
                            bHHO, ...
                            bWOA, ...
                            bGA, ...
                            bBA, ...
                            bPSO, ...
                            bABC, ...
                            bBMWOA, ...
                            bCDLOBA, ...
                            bOBSCA, ...
                            bHHODE, ...
                            bRLHHO, ...
                            bALCPSO, ...
                            bCLPSO, ...
                            bSADE, ...
                        ];
    % select the classifier
    % classifier = 'fknn';
    % classifier = 'svm';
    % classifier = 'kelm';
    classifier = 'knn';
    searchAgentsNum = 20;
    maxFes = 50; % 1000
    runCount = 10;
    folds = 10;

end

%% 原框架内容  需要重新整理——个别存在错误
[row col] = size(algo); %#ok<NCOMMA,ASGLU>
classifiers = {'knn'}; % 分类器 此处可选择的名称有 fknn，svm，kelm，knn，(不完善的地方是，每次分类器只能选择一个)
classifiersNum = size(classifiers, 2);
% classifiersNum1= classifiersNum
%Global parameters
AgentsNum = 20;
%% 迭代次数可进行设置   一般测试可设为：49；正式实验设为：1000
MaxIteration = 50; %1000
NumberOfRuns = 10;
folds = 10; %10
parameters = {'Pop. size', num2cell(AgentsNum), 'MaxIteration', num2cell(MaxIteration), 'NumberOfRuns', num2cell(NumberOfRuns)};
%%%%%%%%% Construct file name %%%%%%%
%% *********** 修改部分 切合栗玉朋_Vesion 1.0 ***************************
kk = 1; % 初始化kk，用于统计框架中需要执行算法的个数

for a = 1:col % 统计可执行算法，并存储在algorithm1中

    if (algo(a) == 1)
        algorithm1(kk) = algorithm(a);
        kk = kk + 1;
    end

end

kk2 = size(algorithm1, 2); % 计算algorithm1中cell数组的个数，用于选择算法以便进行文件命名
algorithm2 = algorithm1{1, kk2}; % 选择最后一个cell数组（注：因为一般会把测试算法放在框架中算法集的最后一个，且该算法一般为新改算法）
data1 = Dataset{1, 1}; % 提取数据集 （注：该版本通过传值的方式，向主函数传递Dataset数据集，且每次只传递一个）
% data2=char(data1);  % 转换为字符串形式
classifiers1 = classifiers{1, 1}; % 提取分类器
% classifiers2=char(classifiers1);
timestr = datestr(now, 'yyyy-mm-dd-HH-MM-SS'); %获取系统时间
timestr1 = datestr(now, 'yyyy-mm-dd');

%% 原文件路径
% dirname = ['result/','',timestr,'-',algorithm{kk}];
%% 改进版本１　（注：测试时可以使用）
% 创建文件存储：dirname = ['result/测试 4.0/','',数据集,'- ',创建时间,'-',目标算法,'-',选择的分类器];
% dirname = ['result/测试 4.0/','',data2,'- ',timestr,'-',algorithm{kk},'-',classifiers2];
%% 最终改进版本　（注：实验版本）
% 创建文件存储：dirname = ['result/实验日期（仅具体到天——通过实验日期对实验结果进行分类）/','',数据集,'- ',创建时间,'-',目标算法,'-',选择的分类器];
dirname = ['result/', timestr1, '-', 'Experiments', '/', '', data1, '- ', timestr, '-', algorithm2, '-', classifiers1];
mkdir(dirname);
%% *********** end of 修改部分 切合栗玉朋_Vesion 1.0 ********************
header = {'Algorithms', 'TransferFun', 'Benchmark', 'Best fitness', 'Error', 'NumOfFeatures', 'Time', 'Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'MCC', 'F-measure'};
header_summary = {'Algorithms', 'TransferFun', 'Benchmark', 'Metric', 'Best fitness', 'Error', 'NumOfFeatures', 'Time'};
header_summary1 = {'Algorithms', 'TransferFun', 'Benchmark', 'Metric', 'Best fitness', 'Error', 'NumOfFeatures', 'Time', 'Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'MCC', 'F-measure'};
header_summary3 = {'Algorithms', 'TransferFun', 'Benchmark', 'Metric', 'Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'MCC', 'F-measure'};
metrics_labels = {'std', 'max', 'min', 'medium'};
x = fix(clock);
str = strtrim(cellstr(num2str(x'))');
strs_spaces = sprintf('-%s', str{:});
trimmed = strtrim(strs_spaces);
% filename=strcat('Experiments-',trimmed);
name = strcat('-Experiments-', trimmed);
filename1 = strcat(classifiers{1}, name);
filename = [dirname, '/', filename1];
filenameALL = strcat(filename, '-AllResults.csv');
filenameConvergences = strcat(filename, '-Convergences.csv');
filenameConvergencesSummary = strcat(filename, '-AvgConvergences.csv');
filenameSummary = strcat(filename, '-Summary.csv');
filenameSummary2 = strcat(filename, '-Summary2.csv');
filenameSummary3 = strcat(filename, '-Summary3.csv');
filenameParameters = strcat(filename, '-Parameters.csv');
filenameBestSolutions = strcat(filename, '-bestSolutions.csv');
filenameFeaturesImportance = strcat(filename, '-FeaturesImportance.csv');
% filenameNumberOfNeurons= strcat(filename,'-NumberOfNeuronesCurve.csv');
filenameFtest = strcat(filename, '-Ftest.csv');
filenameOrder = strcat(filename, '-Order.csv');

dlmcell(filenameALL, header, ',', '-a');
dlmcell(filenameSummary, header_summary1, ',', '-a');
dlmcell(filenameSummary2, header_summary, ',', '-a');
dlmcell(filenameSummary3, header_summary3, ',', '-a');
dlmcell(filenameParameters, parameters, ',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for classifierNum = 1:classifiersNum

    ii = 1; % Counter for total experiments (internal use)

    tic

    for tf = 5
        TFid = tf;
        %Transfer Function Selection  1-4: Sigmoid, 5-8: V-shaped
        %Check the transferFun.m for more details
        for a = 1:col % This loop passes over all selected algorithms

            if (algo(a) == 1)

                for d = 1:size(Dataset, 1) % This loop passes over all datasets (you can make it for d=1:23 to loop over all 23 functions)
                    %                     for jj=1:1

                    %xlRange = 'A1';
                    %xlswrite(filename,header,'Results',xlRange)
                    con = 1;

                    Covergences = [];
                    bestResults = [];
                    bestSolutionss = [];

                    %fn2=[cell2mat(Dataset(d)) '_train'];
                    %A=load([ fn2]);
                    %nVar=size(A,2)-1;
                    %             r=randperm(size(A,1));
                    %             trn=r(1:floor(length(r)/2));
                    %             vald=r(floor(length(r)/2)+1:end);

                    file = strcat(Dataset{d}, '.dat');
                    data = load(file);
                    target = data(:, end);
                    nVar = size(data, 2) - 1;
                    len = size(data, 1);
                    A = data;

                    display(['--------- ', algorithm{a}, ' ---------']);
                    display(['--------- ', Dataset{d}, ' ---------']);
                    %                         fprintf('Run ==> %i \r', i);
                    %------------
                    cvFolds = crossvalind('Kfold', target, folds); %# get indices of 10-fold CV，crossvalind 函数是 MATLAB 中用于生成交叉验证索引的函数。交叉验证是一种用于评估机器学习模型性能的技术，通常用于分割数据集以进行训练和测试。crossvalind 可以帮助你生成交叉验证的索引，以便划分数据集。
                    testIdx = {};

                    for k = 1:folds %floor(len*0.3)   create array of training and testing folds
                        %------------
                        %############ Handle Folds
                        testIdx{k} = (cvFolds == k); %# get indices of test instances

                        %############
                    end

                    %parfor (k=1:folds)
                    parfor (k = 1:folds)
                        display(['Fold: ', num2str(k), ' ---------']);
                        %------------
                        %############ Handle Folds
                        testIx = testIdx{k};
                        trainIdx = ~testIx; %# get indices training instances
                        %                             TrainingData_File=data(trainIdx,:);
                        %                             TestingData_File=data(testIx,:);
                        [findonetest, ~] = find(testIx == 1);
                        [findonetrain, ~] = find(trainIdx == 1);
                        %############
                        % neurons is the number of neurons returened by OP-ELM
                        %                             [TargetFitness,TargetPosition,convergence,acc, Time,cmtest,Cost,Gamma]=optimizeall(a,AgentsNum,MaxIteration,nVar,A,TrainingData_File,TestingData_File,TFid);
                        classifierName = classifiers{classifierNum};
                        [classifierFhd] = Get_Classifiers(classifierName);
                        %% 可以修改算法集：optimizeall  该文件存储的是框架中已有算法的列表
                        % 存储顺序和maincrossvalidationClassifier_v1.文件中algorithm和alg的算法顺序一致，否则可能会出现错位调用，导致实验结果不符合实际
                        % [TargetFitness,TargetPosition,convergence,acc,acc1, sen, spe,pre,mcc,F_measure,Time]=optimizeall2(a,AgentsNum,MaxIteration,nVar,A,findonetrain,findonetest,TFid,classifierFhd);
                        [TargetFitness, TargetPosition, convergence, acc, acc1, sen, spe, pre, mcc, F_measure, Time] = optimizeall_v1(a, AgentsNum, MaxIteration, nVar, A, findonetrain, findonetest, TFid, classifierFhd);
                        %                             ACC(k) = acc1;
                        %                             sens(k) = sen;
                        %                             spec(k) = spe;
                        %                             pres(k) = pre;
                        %                             MCC(k) = mcc;
                        %                             F_means(k) = F;
                        Covergences(k, :) = convergence;
                        %NumberOfNeuronsCurvePrint(k,:)=NumberOfNeuronsCurve;
                        redDim = sum(TargetPosition(:)); %相当于核心集
                        bestResults(k, :) = [TargetFitness acc redDim Time];
                        bestResults1(k, :) = [acc1 sen spe pre mcc F_measure];

                        bestSolutionss(k, :) = TargetPosition;
                        AlgoLabel{k} = algorithm{a};
                        BenchMarkLabel{k} = Dataset{d};
                        TranFunc{k} = TFid;
                    end

                    % save the results of the current run
                    AllResults = [AlgoLabel' TranFunc' BenchMarkLabel' num2cell(bestResults) num2cell(bestResults1)];
                    dlmcell(filenameALL, AllResults, ',', '-a');

                    % save the convergence curves of the current run
                    ConvergenceResults = [AlgoLabel' TranFunc' BenchMarkLabel' num2cell(Covergences)];
                    dlmcell(filenameConvergences, ConvergenceResults, ',', '-a');

                    %Save all best solutions
                    AllBestSolutions = [AlgoLabel' TranFunc' BenchMarkLabel' num2cell(bestSolutionss)];
                    dlmcell(filenameBestSolutions, AllBestSolutions, ',', '-a');

                    % Calculate and save Features importance
                    FeaturesImportance = [AlgoLabel{1}' BenchMarkLabel{1}' num2cell(sum(bestSolutionss, 1) / NumberOfRuns)];
                    dlmcell(filenameFeaturesImportance, FeaturesImportance, ',', '-a');

                    %calculate and save the summary results of this run
                    summarys = [bestResults bestResults1];
                    %                         summary1=mean(bestResults);
                    summary1 = mean(summarys);
                    summaryResults = [AlgoLabel{1}' TranFunc{1}' BenchMarkLabel{1}' 'mean' num2cell(summary1)];
                    dlmcell(filenameSummary, summaryResults, ',', '-a'); % This prints the mean

                    %calculate and save the summary2 results of this run
                    summary2 = [std(bestResults); max(bestResults); min(bestResults); median(bestResults)];
                    %                         summaryResults2=  [cell(4,2) metrics_labels' num2cell(summary2)];
                    summaryResults2 = [AlgoLabel(1:4)' TranFunc(1:4)' BenchMarkLabel(1:4)' metrics_labels' num2cell(summary2)];
                    dlmcell(filenameSummary2, summaryResults2, ',', '-a'); % This prints the mean

                    %calculate and save the summary2 results of this run
                    summary3 = [std(bestResults1); max(bestResults1); min(bestResults1); median(bestResults1)];
                    summaryResults3 = [AlgoLabel(1:4)' TranFunc(1:4)' BenchMarkLabel(1:4)' metrics_labels' num2cell(summary3)];
                    dlmcell(filenameSummary3, summaryResults3, ',', '-a'); % This prints the mean

                    % Calculate Average Convergence Curves for all algorithms andlmcell(filenameSummary2,summaryResults2,',','-a'); % This prints the mean
                    % store them in one file % Ready for plots
                    ConvergenceAvgResults = [AlgoLabel{1}' BenchMarkLabel{1}' num2cell(mean(Covergences, 1))];
                    dlmcell(filenameConvergencesSummary, ConvergenceAvgResults, ',', '-a');

                    %                     end
                end

            end

        end

    end

end

tf = 1; % tf 必须为1 numofrun是每个函数重复运行的次数   d是指有多少个测试集
FridTest_fit1(filenameALL, filenameFtest, tf, d, NumberOfRuns);
FridTest_acc1(filenameALL, filenameFtest, tf, d, NumberOfRuns);
FridTest_fea1(filenameALL, filenameFtest, tf, d, NumberOfRuns);
FridTest_time1(filenameALL, filenameFtest, tf, d, NumberOfRuns);
FridTest_Accuracy(filenameALL, filenameFtest, tf, d, NumberOfRuns);
FridTest_sens(filenameALL, filenameFtest, tf, d, NumberOfRuns);
FridTest_spec(filenameALL, filenameFtest, tf, d, NumberOfRuns);
FridTest_prec(filenameALL, filenameFtest, tf, d, NumberOfRuns);
FridTest_MCC(filenameALL, filenameFtest, tf, d, NumberOfRuns);
FridTest_Fmeasure(filenameALL, filenameFtest, tf, d, NumberOfRuns);
Order_fit1(filenameSummary, filenameOrder, tf, d);
Order_acc1(filenameSummary, filenameOrder, tf, d);
Order_fea1(filenameSummary, filenameOrder, tf, d);
Ordertf_time1(filenameSummary, filenameOrder, tf, d);
Order_Accuracy(filenameSummary, filenameOrder, tf, d);
Order_sens(filenameSummary, filenameOrder, tf, d);
Order_spec(filenameSummary, filenameOrder, tf, d);
Order_prec(filenameSummary, filenameOrder, tf, d);
Order_MCC(filenameSummary, filenameOrder, tf, d);
Order_Fmeasure(filenameSummary, filenameOrder, tf, d);
toc
