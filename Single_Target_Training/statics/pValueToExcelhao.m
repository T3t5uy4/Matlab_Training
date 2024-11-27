function pValueToExcelhao(xlsfilename, fold2)
    %UNTITLED �˴���ʾ�йش˺�����ժҪ
    %   �˴���ʾ��ϸ˵��
    readFilename = xlsfilename;
    writeFilename = xlsfilename;
    [~, ~, rawdata] = xlsread(xlsfilename, 'overall');
    num = size(rawdata);
    flag = 0;
    cur = 1;

    for i = 2:num(1)

        if flag == 0
            cur = cur + 1;
        else
            break;
        end

        if i + 1 <= num(1) && ~strcmp(rawdata(i, 1), rawdata(i + 1, 1))
            flag = 1;
        end

    end

    functioon_num = (num(1) - 1) / (cur - 1);
    algorithm_num = cur - 1;
    fold = fold2;

    sheet_name = 'pValue';
    tag = ['B3:B', num2str(1 + algorithm_num)];
    [~, label] = xlsread(readFilename, 'overall', tag);
    label = label';
    xlswrite(writeFilename, label, sheet_name, 'B1') %����pValue���еĵ�һ��

    for i = 1:functioon_num
        Function_name = ['F', num2str(i)];
        Function_name = {Function_name};
        xlswrite(writeFilename, Function_name, sheet_name, ['A', num2str(i + 1)]); %pValue���е�һ��
    end

    for i = 1:functioon_num
        temp = xlsread(readFilename, ['F', num2str(i)], '', 'basic');
        data(:, i) = temp(:, end);
    end

    pvalue = zeros(functioon_num, (algorithm_num - 1));

    for i = 1:functioon_num

        for j = 1:(algorithm_num - 1)
            pvalue(i, j) = signrank(data(1:fold, i), data(j * fold + 1:j * fold + fold, i));
        end

    end

    xlswrite(writeFilename, pvalue, sheet_name, 'B2') %����F1�е����в������ݣ���C1��ʼ

    % �ϲ�pValue��result
    for i = 1:algorithm_num - 1

        number = 'C' + (i * 3 - 2) - 'A';
        numberOfA = floor(number / 26);
        numberOfVar = mod(number, 26);

        frontName = [];

        if numberOfA ~= 0

            for j = 1:numberOfA
                frontName = [frontName, 'A'];
            end

        end

        behindName = char('A' + numberOfVar);

        xlswrite(writeFilename, pvalue(:, i), 'result & pValue', [[frontName, behindName], '2']);
        xlswrite(writeFilename, {'pvalue'}, 'result & pValue', [[frontName, behindName], '1']);
    end

    % ͳ��ǿ��
    xlswrite(writeFilename, {'CountOfBetter'}, 'result & pValue', ['A', num2str(functioon_num + 4)]);
    xlswrite(writeFilename, {'CountOfWorse'}, 'result & pValue', ['A', num2str(functioon_num + 5)]);
    xlswrite(writeFilename, {'CountOfEqual'}, 'result & pValue', ['A', num2str(functioon_num + 6)]);
    tmpAlpha = ' ';
    tmpResult = repmat(tmpAlpha, functioon_num, algorithm_num - 1);
    tmpdata = xlsread(writeFilename, 'result & pValue');

    for i = 1:algorithm_num - 1
        cntOfBetter = 0;
        cntOfWorse = 0;

        for j = 1:functioon_num

            if tmpdata(j, i * 3) < 0.05

                if tmpdata(j, 1) < tmpdata(j, i * 3 - 1)
                    tmpResult(j, i) = '+';
                    cntOfBetter = cntOfBetter + 1;
                else
                    tmpResult(j, i) = '-';
                    cntOfWorse = cntOfWorse + 1;
                end

            end

        end

        cntOfEqual = functioon_num - cntOfBetter - cntOfWorse;

        number = 'C' + (i * 3 - 1) - 'A';
        numberOfA = floor(number / 26);
        numberOfVar = mod(number, 26);

        frontName = [];

        if numberOfA ~= 0

            for j = 1:numberOfA
                frontName = [frontName, 'A'];
            end

        end

        behindName = char('A' + numberOfVar);

        xlswrite(writeFilename, tmpResult(:, i), 'result & pValue', [[frontName, behindName], '2']);
        xlswrite(writeFilename, cntOfBetter, 'result & pValue', [[frontName, behindName], num2str(functioon_num + 4)]);
        xlswrite(writeFilename, cntOfWorse, 'result & pValue', [[frontName, behindName], num2str(functioon_num + 5)]);
        xlswrite(writeFilename, cntOfEqual, 'result & pValue', [[frontName, behindName], num2str(functioon_num + 6)]);
    end
