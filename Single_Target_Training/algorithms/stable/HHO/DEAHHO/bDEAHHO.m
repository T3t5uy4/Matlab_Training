% The Harris hawks optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve, time] = bDEAHHO(searchAgentsNum, maxFes, A, trn, vald, TFid, classifierFhd)
    % Initialize position vector and fitness for the best
    tic
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, 1, 0);
    fitness = inf * ones(searchAgentsNum, 1);
    convergenceCurve = [];

    t = 0;
    F = 0.5; % Scaling factor
    ub = 1;
    lb = 0;

    while t < maxFes

        for i = 1:size(positions, 1)
            % Fitness of locations
            fitness(i) = AccSz2(positions(i, :), A, trn, vald, classifierFhd);

            if fitness(i) < bestFitness
                bestFitness = fitness(i);
                bestPosition = positions(i, :);
            end

        end

        % Update the Harris hawk's position
        for i = 1:size(positions, 1)
            % Get the escape energies
            E0 = 2 * rand - 1;
            E = 2 * E0 * (1 - fe / maxFes);

            if abs(E) >= 1 % Exploration stage
                q = rand;
                randHawkIdx = floor(searchAgentsNum * rand + 1);
                randPosition = positions(randHawkIdx, :);

                if q < 0.5
                    % Perch based on other family members
                    positions(i, :) = randPosition - rand * abs(randPosition - 2 * rand * positions(i, :));

                    for i = 1:dim
                        positions(i, j) = transferFun(positions(i, j), positions(i, j), TFid);
                    end

                elseif q >= 0.5
                    % Perch on a random tall tree (random site inside group's home range)
                    positions(i, :) = bestPosition - mean(positions) - rand * ((ub - lb) * rand + lb);

                    for j = 1:dim
                        positions(i, j) = transferFun(positions(i, j), positions(i, j), TFid);
                    end

                end

            elseif abs(E) < 1 % Development stage
                r = rand;

                if r >= 0.5 && abs(E) < 0.5 % Hard besiege
                    positions(i, :) = bestPosition - E * abs(bestPosition - positions(i, :));

                    for j = 1:dim
                        positions(i, j) = transferFun(positions(i, j), positions(i, j), TFid);
                    end

                elseif r >= 0.5 && abs(E) >= 0.5 % Soft besiege
                    jumpStrength = 2 * (1 - rand); % Random jump strength
                    positions(i, :) = bestPosition - positions(i, :) - E * abs(jumpStrength * bestPosition - positions(i, :));

                    for j = 1:dim
                        positions(i, j) = transferFun(positions(i, j), positions(i, j), TFid);
                    end

                elseif r < 0.5 && abs(E) >= 0.5 % Soft besiege with motions
                    jumpStrength = 2 * (1 - rand); % Random jump strength
                    position1 = bestPosition - E * abs(jumpStrength * bestPosition - positions(i, :));

                    for j = 1:dim
                        position1(1, j) = transferFun(position1(1, j), position1(1, j), TFid);
                    end

                    fitness1 = AccSz2(position1, A, trn, vald, classifierFhd);

                    if fitness1 < fitness(i)
                        fitness(i) = fitness1;
                        positions(i, :) = position1;
                    else
                        position2 = bestPosition - E * abs(jumpStrength * bestPosition - positions(i, :)) + rand(1, dim) .* Levy(dim);

                        for j = 1:dim
                            position2(1, j) = transferFun(position2(1, j), position2(1, j), TFid);
                        end

                        fitness2 = AccSz2(position2, A, trn, vald, classifierFhd);

                        if fitness2 < fitness(i)
                            fitness(i) = fitness2;
                            positions(i, :) = position2;
                        end

                    end

                elseif r < 0.5 && abs(E) < 0.5 % Hard besiege with motions
                    jumpStrength = 2 * (1 - rand); % Random jump strength
                    position1 = bestPosition - E * abs(jumpStrength * bestPosition - mean(positions));

                    for j = 1:dim
                        position1(1, j) = transferFun(position1(1, j), position1(1, j), TFid);
                    end

                    fitness1 = AccSz2(position1, A, trn, vald, classifierFhd);

                    if fitness1 < fitness(i)
                        fitness(i) = fitness1;
                        positions(i, :) = position1;
                    else
                        position2 = bestPosition - E * abs(jumpStrength * bestPosition - mean(positions)) + rand(1, dim) .* Levy(dim);

                        for j = 1:dim
                            position2(1, j) = transferFun(position2(1, j), position2(1, j), TFid);
                        end

                        fitness2 = AccSz2(position2, A, trn, vald, classifierFhd);

                        if fitness2 < fitness(i)
                            fitness(i) = fitness2;
                            positions(i, :) = position2;
                        end

                    end

                end

            end

            r = rand;
            feNorm = fe / maxFes;
            p1 = (exp(-feNorm) - exp(feNorm) + exp(1) - exp(-1)) / 3;
            p2 = 1 - (exp(- (feNorm - 0.5)) + exp(feNorm - 0.5)) / 3;
            CR = 1 - 0.3 * (1 - fe / maxFes);

            if r <= p1
                % DE_rand
                % Select three random indices
                indices = randperm(searchAgentsNum, 3);
                positionRand1 = positions(indices(1), :);
                positionRand2 = positions(indices(2), :);
                positionRand3 = positions(indices(3), :);

                % Mutation
                mutant = positionRand1 + (positionRand2 - positionRand3) .* Levy(dim);

                for j = 1:dim
                    mutant(1, j) = mutant(positions(1, j), mutant(1, j), TFid);
                end

            elseif r <= p1 + p2
                % DE_currentToBest
                % Select two random indices
                indices = randperm(searchAgentsNum, 2);
                positionRand1 = positions(indices(1), :);
                positionRand2 = positions(indices(2), :);

                % Mutation
                mutant = positions(i, :) + (positionRand1 - positionRand2) .* Levy(dim);

                for j = 1:dim
                    mutant(1, j) = mutant(positions(1, j), mutant(1, j), TFid);
                end

            else
                % DE_currentToBest
                % Select two random indices
                indices = randperm(searchAgentsNum, 2);
                positionRand1 = positions(indices(1), :);
                positionRand2 = positions(indices(2), :);

                % Mutation
                mutant = positions(i, :) + F * (bestPosition - positions(i, :)) + (positionRand1 - positionRand2) .* Levy(dim);

                for j = 1:dim
                    mutant(1, j) = mutant(positions(1, j), mutant(1, j), TFid);
                end

            end

            % Crossover
            trial = positions(i, :);
            j_rand = randi(dim); % Random index for crossover

            for j = 1:dim

                if rand < CR || j == j_rand
                    trial(j) = mutant(j);
                end

            end

            % Selection
            trialFitness = trnasferFun(trial(1, j), trial(1, j), TFid);

            if trialFitness < fitness(i)
                positions(i, :) = trial;
                fitness(i) = trialFitness;
            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

    Time = toc;
end
