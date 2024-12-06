% last modify by luojie in 20171018
function [lb, Positionbest, Convergence_curve] = GA(SearchAgents_no, MaxFEs, lb, ub, dim, fobj)
    %% ������ʼ��

    ga_option = struct('maxgen', MaxFEs, 'sizepop', SearchAgents_no, 'pCrossover', 0.3, 'pMutation', 0.1, ...
        'cbound', [lb, ub], 'gbound', [lb, ub], 'v', 3);

    FEs = 0;
    c_len_chromosome = ceil(log2((ga_option.cbound(2) - ga_option.cbound(1)) * 100));
    g_len_chromosome = ceil(log2((ga_option.gbound(2) - ga_option.gbound(1)) * 100));
    len_chromosome = c_len_chromosome + g_len_chromosome;

    % ����Ⱥ��Ϣ����Ϊһ���ṹ��
    individuals = struct('fitness', zeros(1, ga_option.sizepop), ...
        'chromosome', zeros(ga_option.sizepop, len_chromosome));
    % ÿһ����Ⱥ��ƽ����Ӧ��
    avgfitness_gen = zeros(1, ga_option.maxgen);
    % ÿһ����Ⱥ�������Ӧ��
    bestfitness_gen = [];
    Convergence_curve = [];

    %% ��ʼ����Ⱥ
    fobjvalue = zeros(1, ga_option.sizepop);

    for i = 1:ga_option.sizepop
        % ����
        individuals.chromosome(i, :) = unidrnd(2, 1, len_chromosome) - 1;
        % ����
        [c, g] = ga_decode(individuals.chromosome(i, :), ga_option.cbound, ga_option.gbound);
        % �����ʼ��Ӧ��(CV׼ȷ��)
        fobjvalue(i) = fobj([c, g]);
        FEs = FEs + 1;
    end

    globleMax = max(fobjvalue);
    globleMin = min(fobjvalue);
    individuals.fitness = getFitness(fobjvalue, globleMax, globleMin); %������Ĵ�����p
    % ����ѵ���Ӧ�Ⱥ���õ�Ⱦɫ���λ��
    [bestfitness, bestindex] = max(individuals.fitness);
    % ��õ�Ⱦɫ��
    bestchromosome = individuals.chromosome(bestindex, :);
    t = 1;
    %% ����Ѱ��
    while FEs < ga_option.maxgen %i=1:ga_option.maxgen
        % Selection Operator
        individuals = Selection(individuals, ga_option);
        % Crossover Operator
        individuals = Crossover(individuals, ga_option);
        % Mutation Operator
        individuals = Mutation(individuals, ga_option);

        % ������Ӧ��
        fobjvalue = zeros(1, ga_option.sizepop);

        for j = 1:ga_option.sizepop
            % ����
            [c, g] = ga_decode(individuals.chromosome(j, :), ga_option.cbound, ga_option.gbound);
            fobjvalue(j) = fobj([c, g]);
            FEs = FEs + 1;
        end

        localMax = max(fobjvalue);
        localMin = min(fobjvalue);

        if localMax > globleMax,
            globleMax = localMax;
        end

        if localMin < globleMin,
            globleMin = localMin;
        end

        individuals.fitness = getFitness(fobjvalue, globleMax, globleMin); %������Ĵ�����p
        % ����ѵ���Ӧ�Ⱥ���õ�Ⱦɫ���λ��
        [new_bestfitness, bestindex] = max(individuals.fitness);
        % ��õ�Ⱦɫ��
        new_bestchromosome = individuals.chromosome(bestindex, :);

        if new_bestfitness >= bestfitness
            bestfitness = new_bestfitness;
            bestchromosome = new_bestchromosome;
        end

        % ��һ��Ⱦɫ��������Ӧ��
        bestfitness_gen(t) = bestfitness;
        % ��һ��Ⱦɫ���ƽ����Ӧ��
        avgfitness_gen(t) = sum(individuals.fitness) / ga_option.sizepop;
        %for output
        Positionbest = ga_decode(bestchromosome, ga_option.cbound, ga_option.gbound);
        bestValue = fobj(Positionbest);
        FEs = FEs + 1;

        if t > 1 && bestValue > Convergence_curve(t - 1)
            Convergence_curve(t) = Convergence_curve(t - 1);
        else
            Convergence_curve(t) = bestValue;
        end

        t = t + 1;
    end

end

%% sub function ga_decode
function [c, g] = ga_decode(chromosome, cbound, gbound)
    % ga_decode by faruto
    % Email:farutoliyang@gmail.com
    % 2009.10.08
    c_len_chromosome = ceil(log2((cbound(2) - cbound(1)) * 100));
    g_len_chromosome = ceil(log2((gbound(2) - gbound(1)) * 100));
    len_chromosome = c_len_chromosome + g_len_chromosome;

    cdec = bin2dec(num2str(chromosome(1:c_len_chromosome)));
    gdec = bin2dec(num2str(chromosome(c_len_chromosome + 1:len_chromosome)));

    c = cbound(1) + cdec * (cbound(2) - cbound(1)) / (2 ^ (c_len_chromosome) - 1);
    g = gbound(1) + gdec * (gbound(2) - gbound(1)) / (2 ^ (g_len_chromosome) - 1);
end

%% sub function Selection
function individuals_afterSelect = Selection(individuals, ga_option)
    % Selection by faruto
    % Email:farutoliyang@gmail.com
    % 2009.10.08
    individuals_afterSelect = individuals;
    sum_fitness = sum(individuals.fitness);
    P = individuals.fitness / sum_fitness;
    Q = zeros(1, ga_option.sizepop);

    for k = 1:ga_option.sizepop
        Q(k) = sum(P(1:k));
    end

    for i = 1:ga_option.sizepop
        r = rand;

        while r == 0
            r = rand;
        end

        k = 1;

        while k <= ga_option.sizepop - 1 && r > Q(k)
            k = k + 1;
        end

        %     individuals_afterSelect.fitness(i) = individuals.fitness(k);
        individuals_afterSelect.chromosome(i, :) = individuals.chromosome(k, :);
    end

end

%% sub function Crossover
function individuals_afterCross = Crossover(individuals, ga_option)
    % Crossover by faruto
    % Email:farutoliyang@gmail.com
    % 2009.10.08
    individuals_afterCross = individuals;
    c_len_chromosome = ceil(log2((ga_option.cbound(2) - ga_option.cbound(1)) * 100));
    g_len_chromosome = ceil(log2((ga_option.gbound(2) - ga_option.gbound(1)) * 100));
    len_chromosome = c_len_chromosome + g_len_chromosome;

    for i = 1:ga_option.sizepop
        % ������ʾ����Ƿ���н���
        r = rand;

        if r > ga_option.pCrossover
            continue;
        end

        % ���ѡ������Ⱦɫ����н���
        pick = rand(1, 2);

        while prod(pick) == 0
            pick = rand(1, 2);
        end

        index = ceil(pick .* ga_option.sizepop);

        % ���ѡ�񽻲�λ��
        pos_cross = unidrnd(len_chromosome - 1);
        % ���н���
        individuals_afterCross.chromosome(index(1), pos_cross + 1:len_chromosome) ...
            = individuals.chromosome(index(2), pos_cross + 1:len_chromosome);
        individuals_afterCross.chromosome(index(2), pos_cross + 1:len_chromosome) ...
            = individuals.chromosome(index(1), pos_cross + 1:len_chromosome);
    end

end

%% sub function Mutation
function individuals_afterMutate = Mutation(individuals, ga_option)
    % Mutation by faruto
    % Email:farutoliyang@gmail.com
    % 2009.10.08
    individuals_afterMutate = individuals;
    c_len_chromosome = ceil(log2((ga_option.cbound(2) - ga_option.cbound(1)) * 100));
    g_len_chromosome = ceil(log2((ga_option.gbound(2) - ga_option.gbound(1)) * 100));
    len_chromosome = c_len_chromosome + g_len_chromosome;

    for i = 1:ga_option.sizepop
        % ������ʾ����Ƿ���н���
        r = rand;

        if r > ga_option.pMutation
            continue;
        end

        % ���ѡ��һ��Ⱦɫ����б���
        pick = unidrnd(ga_option.sizepop);
        % ���ѡ�����λ��
        pos_mutate = unidrnd(len_chromosome);
        % ���б���
        if individuals_afterMutate.chromosome(pick, pos_mutate) == 0
            individuals_afterMutate.chromosome(pick, pos_mutate) = 1;
        else
            individuals_afterMutate.chromosome(pick, pos_mutate) = 0;
        end

    end

end

function fitness = getFitness(objValues, globleMax, globleMin)

    if globleMax == globleMin
        fitness = 0.5 * ones(1, size(objValues, 2));
    else
        fitness = 1 - 0.9 * (objValues - globleMin) / (globleMax - globleMin);
    end

end
