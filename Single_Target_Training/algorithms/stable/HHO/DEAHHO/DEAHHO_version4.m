% The Harris hawks optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = DEAHHO_version4(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    fitness = [];
    convergenceCurve = [];

    fe = 0;
    t = 0;
    F = 0.5;
    CR = 0.9;

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
                elseif q >= 0.5
                    positions(i, :) = bestPosition - mean(positions) - rand * ((ub - lb) * rand + lb);
                    % Perch on a random tall tree (random site inside group's home range)
                end

            elseif abs(E) < 1 % Development stage
                r = rand;

                if r >= 0.5 && abs(E) < 0.5 % Hard besiege
                    positions(i, :) = bestPosition - E * abs(bestPosition - positions(i, :));
                elseif r >= 0.5 && abs(E) >= 0.5 % Soft besiege
                    jumpStrength = 2 * (1 - rand); % Random jump strength
                    positions(i, :) = bestPosition - positions(i, :) - E * abs(jumpStrength * bestPosition - positions(i, :));
                elseif r < 0.5 && abs(E) >= 0.5 % Soft besiege with motions
                    jumpStrength = 2 * (1 - rand); % Random jump strength
                    position1 = bestPosition - E * abs(jumpStrength * bestPosition - positions(i, :));

                    if fobj(position1) < fobj(positions(i, :))
                        positions(i, :) = position1;
                    else
                        position2 = bestPosition - E * abs(jumpStrength * bestPosition - positions(i, :)) + rand(1, dim) .* Levy(dim);

                        if fobj(position2) < fobj(positions(i, :))
                            positions(i, :) = position2;
                        end

                    end

                    fe = fe + 1;

                elseif r < 0.5 && abs(E) < 0.5 % Hard besiege with motions
                    jumpStrength = 2 * (1 - rand); % Random jump strength
                    position1 = bestPosition - E * abs(jumpStrength * bestPosition - mean(positions));

                    if fobj(position1) < fobj(positions(i, :))
                        positions(i, :) = position1;
                    else
                        position2 = bestPosition - E * abs(jumpStrength * bestPosition - mean(positions)) + rand(1, dim) .* Levy(dim);

                        if fobj(position2) < fobj(positions(i, :))
                            positions(i, :) = position2;
                        end

                    end

                    fe = fe + 1;

                end

            end

            p = floor(rand * 8) + 1;

            switch p
                case 1

                    % Select three random indices
                    indices = randperm(searchAgentsNum, 2);
                    positionRand1 = positions(indices(1), :);
                    positionRand2 = positions(indices(2), :);

                    % Mutation
                    mutant = bestPosition + F * (positionRand1 - positionRand2);
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

                case 2

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

                case 3

                    % Select three random indices
                    indices = randperm(searchAgentsNum, 2);
                    positionRand1 = positions(indices(1), :);
                    positionRand2 = positions(indices(2), :);

                    % Mutation
                    mutant = positions(i, :) + F * (positionRand1 - positionRand2);
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

                case 4

                    % Select three random indices
                    indices = randperm(searchAgentsNum, 4);
                    positionRand1 = positions(indices(1), :);
                    positionRand2 = positions(indices(2), :);
                    positionRand3 = positions(indices(3), :);
                    positionRand4 = positions(indices(4), :);

                    % Mutation
                    mutant = positions(i, :) + F * (positionRand1 - positionRand2) + F * (positionRand3 - positionRand4);
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

                case 5

                    % Select three random indices
                    indices = randperm(searchAgentsNum, 2);
                    positionRand1 = positions(indices(1), :);
                    positionRand2 = positions(indices(2), :);

                    % Mutation
                    mutant = positions(i, :) + F * (bestPosition - positions(i, :)) + F * (positionRand1 - positionRand2);
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

                case 6

                    % Select three random indices
                    indices = randperm(searchAgentsNum, 4);
                    positionRand1 = positions(indices(1), :);
                    positionRand2 = positions(indices(2), :);
                    positionRand3 = positions(indices(3), :);
                    positionRand4 = positions(indices(4), :);

                    % Mutation
                    mutant = positions(i, :) + F * (bestPosition - positions(i, :)) + F * (positionRand1 - positionRand2) + F * (positionRand3 - positionRand4);
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

                case 7

                    % Select three random indices
                    indices = randperm(searchAgentsNum, 3);
                    positionRand1 = positions(indices(1), :);
                    positionRand2 = positions(indices(2), :);
                    positionRand3 = positions(indices(3), :);

                    % Mutation
                    mutant = positionRand1 + F * (positionRand2 - positionRand3);
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

                case 8

                    % Select three random indices
                    indices = randperm(searchAgentsNum, 5);
                    positionRand1 = positions(indices(1), :);
                    positionRand2 = positions(indices(2), :);
                    positionRand3 = positions(indices(3), :);
                    positionRand4 = positions(indices(4), :);
                    positionRand5 = positions(indices(5), :);

                    % Mutation
                    mutant = positionRand1 + F * (positionRand2 - positionRand3) + F * (positionRand4 - positionRand5);
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

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end
