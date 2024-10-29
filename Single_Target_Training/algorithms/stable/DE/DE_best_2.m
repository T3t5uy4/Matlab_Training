function [bestFitness, bestPosition, convergenceCurve] = DE_best_2(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    convergenceCurve = [];
    fitness = [];

    F = 0.5; % Scaling factor
    CR = 0.9; % Crossover probability
    fe = 0;
    t = 0;

    while fe < maxFes

        for i = 1:size(positions, 1)
            % Check boundries
            FU = positions(i, :) > ub;
            FL = positions(i, :) < lb;
            positions(i, :) = (positions(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
            % Fitness of locations
            fitness(i) = fobj(positions(i, :));
            fe = fe + 1;

            if fitness(i) < bestFitness
                bestFitness = fitness(i);
                bestPosition = positions(i, :);
            end

        end

        for i = 1:size(positions, 1)
            % Select three random indices
            indices = randperm(searchAgentsNum, 4);
            positionRand1 = positions(indices(1), :);
            positionRand2 = positions(indices(2), :);
            positionRand3 = positions(indices(3), :);
            positionRand4 = positions(indices(4), :);

            % Mutation
            mutant = bestPosition + F * (positionRand1 - positionRand2) + F * (positionRand3 - positionRand4);
            mutant = max(min(mutant, ub), lb); % Ensure within bounds

            % Crossover
            trial = positions(i, :);
            j_rand = randi(dim); % Random index for crossover

            for j = 1:dim

                if rand < CR || j == j_rand
                    trial(j) = mutant(j);
                end

            end

            % Selection
            trialFitness = fobj(trial);
            fe = fe + 1;

            if trialFitness < fitness(i)
                positions(i, :) = trial;
                fitness(i) = trialFitness;
                fe = fe + 1;
            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end
