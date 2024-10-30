%**************************************************************************************************
%  CLPSO�� Comprehensive Learning Particle Swarm Optimizer for Global Optimization of Multimodal Functions
%  We obtained the MATLAB source code from the authors
%  Date: 2014/11/14
%**************************************************************************************************

function [D, Leader_pos, Convergence_curve] = CLPSO(SearchAgents_no, MaxFEs, lb, ub, dim, fobj)

    if size(ub, 2) == 1 %���ؾ���ub������Ϊ1��
        ub = ones(1, dim) * ub;
        lb = ones(1, dim) * lb;
    end

    lu = [lb; ub];

    D = dim;
    Xmin = lu(1, :); %x����С��Χ
    Xmax = lu(2, :); %x�����Χ

    rand('seed', sum(100 * clock));

    % parameters setting for CLPSO
    popsize = SearchAgents_no; %���Ӹ���
    maxGEN = MaxFEs / popsize; %����Ȩ�ع�ʽ��һ����
    iwt = 0.9 - (1:maxGEN) * (0.7 / maxGEN); %����Ȩ��
    c = 1.49445;
    FES = 0;
    % Initialize the main population
    X = repmat(Xmin, popsize, 1) + rand(popsize, D) .* (repmat(Xmax - Xmin, popsize, 1)); %��ʼ���������λ��
    val_X = zeros(size(X, 1), 1); %����D��һ�е�0����

    for i = 1:size(X, 1) %����ÿһ�����ӵ���Ӧֵ
        val_X(i) = fobj(X(i, :));
        FES = FES + 1;
    end

    pBest = X; val_pBest = val_X; %����������ʷ��Сֵ��λ�ú͸����ӵ���Ӧֵ
    [~, indexG] = min(val_pBest); %����Ⱥ��ʷ��Сֵ������λ��
    gBest = pBest(indexG, :); val_gBest = val_pBest(indexG, :); %����Ⱥ��ʷ��Сֵ��λ�ú���С��Ӧֵ
    Vmax = (Xmax - Xmin) * 0.2; Vmin = -Vmax; %��������ٶȺ���С�ٶ�
    V = repmat(Vmin, popsize, 1) + rand(popsize, D) .* repmat(Vmax - Vmin, popsize, 1); %���һ����ʼ�ٶ�

    % Learning Probability Pc
    t = 0:1 / (popsize - 1):1;
    t = 5 .* t;
    Pc = 0.0 + (0.5 - 0.0) .* (exp(t) - exp(t(1))) ./ (exp(t(popsize)) - exp(t(1))); %����ÿ�����ӵ�pcֵ

    % Refreshing Gap
    gapm = 5; stay_num = zeros(1, popsize); %gapm�൱���㷨��m��ֵ��stay_num��¼ÿ�����ӵļ������

    for i = 1:popsize
        pBest_ind(i, :) = LearnIndex_CLPSO(val_pBest, popsize, D, i, Pc(i)); %��ÿ�����ӵ�ά�Ƚ��бȽϣ�ѧϰ����������������ӵľ���
    end

    l = 1;
    GEN = 1; %���������ļ�����
    Convergence_curve = []; % record the best results
    %l=1;
    while FES < MaxFEs

        for i = 1:popsize

            % update exemplar index
            if stay_num(i) > gapm
                pBest_ind(i, :) = LearnIndex_CLPSO(val_pBest, popsize, D, i, Pc(i)); %�ܵ�������������Σ������ӵ���ʷ����ֵ��û�и��£��������ѧϰ��������
                stay_num(i) = 0;
            end

            %  Compehensive Learning Strategy
            for j = 1:D
                pBest_f(i, j) = pBest(pBest_ind(i, j), j); %����������ά���ϵ�ֵ��ֲ����ǰ���ӵ�ά����
            end

            V(i, :) = iwt(GEN) * V(i, :) + c * rand(1, D) .* (pBest_f(i, :) - X(i, :)); % update velocity
            V(i, :) = boundConstraint_absorb(V(i, :), Vmin, Vmax);
            X(i, :) = X(i, :) + V(i, :); % update position

            if all(X(i, :) <= Xmax) && all(X(i, :) >= Xmin) % X(i,:) is feasible
                val = fobj(X(i, :));
                FES = FES + 1;

                if val < val_pBest(i) % update pBest
                    pBest(i, :) = X(i, :); val_pBest(i) = val;
                    stay_num(i) = 0; %������ʷ����ֵ���º󣬼�������Ϊ0

                    if val < val_gBest % update gBestȫ������
                        gBest = X(i, :); val_gBest = val;
                    end

                else
                    stay_num(i) = stay_num(i) + 1; %������ʷ����ֵû�и���ʱ��������һ
                end

            end

        end

        %     convergence = [convergence val_gBest];
        Convergence_curve(l) = val_gBest; %��¼���ε�����ȫ������ֵ
        l = l + 1; %�ܵ���������һ
        GEN = GEN + 1;

        if (GEN >= maxGEN) && (l < MaxFEs)
            GEN = GEN - 1;
        end

    end

    % bestScore=convergence(end);
    Leader_pos = gBest;

end
