clc;
clear;
close all;

% select Dataset Collections
datasets = {};

for i = 1:size(datasets, 2)
    dataset = datasets(i);
    % select the count of the train
    count = 1;

    for j = 1:count
        % select single train or multi train
        single_train(dataset);
        multi_train(dataset);
    end

end

% Dataset = {
%            'Cleveland_heart' %303*14
%            'cmc' % 1473 * 10
%            'CNS' % 60 * 7130 内存不足
%            'CTG3' % 2126 * 22
%            'Dermatology' % 358 * 35
%            'DLBCL' % 77*5470 内存不足
%            'German' % 1000 * 25
%            'glass' % 214 * 10
%            'hepatitisfulldata' % 155 * 20
%            'JPNdata' % 152 * 11
%            'Leukemia1' %72 * 5328
%            'Leukemia2' %72*11226
%            'Lung_Cancer' %203 * 12601
% % 'lungcancer_3class'  % 32 * 57
%            'Prostate_Tumor' %102*10510
%            'SRBCT' %83*2309
%            'thyroid_2class' % 187 * 9
%            'Tumors_9' % 60*5727
%            'Tumors_11' %174*12534
%            'Tumors_14' % 308*15010
%            'USAdata' % 2336 * 11
%            'vehicle' %846*19
%            'wdbc' %569*31
%            'Wielaw' %240*31   95个
%            'blood' %51*23
%            'Brain_Tumor2' %50*10368
%            'JPNdata' % 152 * 11
%            'thyroid' %215*6
%            'Australian' %690*15
%            'CongressEW' % 435 *17
%            'SpectEW' % 267 *23
%            'HeartEW' %270*14
%            'IonosphereEW' %  351*35
%            'Dermatology' % 358 * 35
%            'WineEW' %178*14
%            'transfusion' %748*5
%            'segment' %2310*18
%            'German' %1000*25
%            'SRBCT' %83*2309
%            'Leukemia' %72*7131
%            'KrvskpEW' %3196*37
%            'Lung_Cancer' %203 * 12601
%            'teachingassistant' %151*6
%            'thyroid_2class' %187*9
%            'vehicle' %846*19
%            'glass' % 214 * 10
%            'Cleveland_heart' %303*14
%            'Vote' % 101*17
%            'Wielaw' %240*31
%            'BreastEW' % 569*31
%            'DLBCL' %77*5400
%            'heart' % 270*14
%            'penglungEW' % 73*326
%            'iris' %150*5
%            'clean1' %476*167
%            'Prostate_Tumor' %102*10510
%            'WaveformEW' %  5000*41
%            'semeion' %1593*266
%            'USAdata' %2336*11
%            'zoo' %101*17
%            'Breastcancer' % 699*10
%            'Lymphography' % 148*19
%            'primary-tumor' %339*18
%            'Tic-tac-toe' %958*10
%            'Breastcancer' % 699*10
%            'BreastEW' % 569*31
% %'Exactly' % 1000*14
% % 'Exactly2'  % 1000*14
%            'HeartEW' %270*14
%            'primary-tumor' % 339*18
%            'heart' % 270*14
%            'Parkinson' %195*23
%            'IonosphereEW' %  351*35
%            'Lymphography' % 148*19
%            'M-of-n' % 100*14
%            'penglungEW' % 73*326
%            'SonarEW' % 208*61
%            'SpectEW' % 267 *23
%            'CongressEW' % 435 *17
%            'KrvskpEW' % 208*61
%            'Tic-tac-toe' %958*10
%            'Vote' % 101*17
%            'WaveformEW' %  5000*41
%            'Zoo' %101*17
%            'clean1' %476*167
%            'clean2' %6598*167
%            'semeion' %1593*266
% % 'Colon'   %62*high dimension   62*2001
%            'Leukemia' %72*7131  high dimension 内存不足
%            'Australian' %690*15
%            'Brain_Tumor1' %90*5921 内存不足
%            'Brain_Tumor2' %50*10368  内存不足
%            };

% %% 选择数据集——控制起始位置以及选择特定的数据集
% k = 1; % 数据集起始位置
% % 测试参数
% K = 2;
% % 实验参数
% % K=size(Dataset,１); % 获取数据集集合中数据集的个数，用于控制数据集的选择
% % K=22; % 数据集的总个数：可以在已有数据集中自定义
% %%
% while k <= K %（通过选择多个数据集经行单独实验，来控制实验次数）　　ｊ表示主函数被执行的次数　
%     %% 控制实验次数 1
%     i = 0;
%     % 测试参数　　
%     j = 1;
%     % 实验参数
%     % j = 10;
%     % 选择数据集
%     dataset = Dataset(k); % 依次选择数据集进行实验
%     %% 控制每个数据集在实验时执行的次数
%     while i < j % 控制执行次数 2  空子的是每个数据集被执行实验的次数
%         maincrossvalidationClassifier_v1(dataset); % 调取实验主函数，执行实验
%         i = i + 1; % 执行次数控制
%     end

%     k = k + 1; % 用于控制选择数据集
% end
