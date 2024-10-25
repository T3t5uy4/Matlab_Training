function [bestf, bestx, BestCost] = DE(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction)
    %% DE Parameters

    % nVar= Number of Decision Variables
    VarSize = [1 nVar]; % Decision Variables Matrix Size
    % VarMin=Lower Bound of Decision Variables
    % VarMax= Upper Bound of Decision Variables
    % MaxIt= Maximum Number of Iterations
    % nPop= Population Size
    beta_min = 0.2; % Lower Bound of Scaling Factor
    beta_max = 0.8; % Upper Bound of Scaling Factor
    pCR = 0.2; % Crossover Probability
    %% Initialization
    empty_individual.Position = [];
    empty_individual.Cost = [];
    BestSol.Cost = inf;
    pop = repmat(empty_individual, nPop, 1);

    for i = 1:nPop
        pop(i).Position = unifrnd(VarMin, VarMax, VarSize);

        pop(i).Cost = CostFunction(pop(i).Position);

        if pop(i).Cost < BestSol.Cost
            BestSol = pop(i);
        end

    end

    BestCost = zeros(MaxIt, 1);
    %% DE Main Loop
    for it = 1:MaxIt

        for i = 1:nPop

            x = pop(i).Position;

            A = randperm(nPop);

            A(A == i) = [];

            a = A(1);
            b = A(2);
            c = A(3);

            % Mutation
            %beta=unifrnd(beta_min,beta_max);
            beta = unifrnd(beta_min, beta_max, VarSize);
            y = pop(a).Position + beta .* (pop(b).Position - pop(c).Position);
            y = max(y, VarMin);
            y = min(y, VarMax);

            % Crossover
            z = zeros(size(x));
            j0 = randi([1 numel(x)]);

            for j = 1:numel(x)

                if j == j0 || rand <= pCR
                    z(j) = y(j);
                else
                    z(j) = x(j);
                end

            end

            NewSol.Position = z;
            NewSol.Cost = CostFunction(NewSol.Position);

            if NewSol.Cost < pop(i).Cost
                pop(i) = NewSol;

                if pop(i).Cost < BestSol.Cost
                    BestSol = pop(i);
                end

            end

        end

        % Update Best Cost
        BestCost(it) = BestSol.Cost;

    end

    %% Show Results

    bestx = BestSol.Position;
    bestf = BestCost(end);
end
