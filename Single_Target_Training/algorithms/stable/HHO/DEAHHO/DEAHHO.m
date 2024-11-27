% The Harris hawks optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = DEAHHO(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    fitness = [];
    convergenceCurve = [];

    t = 0;
    fe = 0;
    F = 0.5; % Scaling factor

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
                    % Perch on a random tall tree (random site inside group's home range)
                    positions(i, :) = bestPosition - mean(positions) - rand * ((ub - lb) * rand + lb);
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
                    fitness1 = fobj(position1);
                    fe = fe + 1;

                    if fitness1 < fitness(i)
                        fitness(i) = fitness1;
                        positions(i, :) = position1;
                    else
                        position2 = bestPosition - E * abs(jumpStrength * bestPosition - positions(i, :)) + rand(1, dim) .* Levy(dim);
                        fitness2 = fobj(position2);
                        fe = fe + 1;

                        if fitness2 < fitness(i)
                            fitness(i) = fitness2;
                            positions(i, :) = position2;
                        end

                    end

                elseif r < 0.5 && abs(E) < 0.5 % Hard besiege with motions
                    jumpStrength = 2 * (1 - rand); % Random jump strength
                    position1 = bestPosition - E * abs(jumpStrength * bestPosition - mean(positions));
                    fitness1 = fobj(position1);
                    fe = fe + 1;

                    if fitness1 < fitness(i)
                        fitness(i) = fitness1;
                        positions(i, :) = position1;
                    else
                        position2 = bestPosition - E * abs(jumpStrength * bestPosition - mean(positions)) + rand(1, dim) .* Levy(dim);
                        fitness2 = fobj(position2);
                        fe = fe + 1;

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
                mutant = max(min(mutant, ub), lb); % Ensure within bounds
            elseif r <= p1 + p2
                % DE_currentToBest
                % Select two random indices
                indices = randperm(searchAgentsNum, 2);
                positionRand1 = positions(indices(1), :);
                positionRand2 = positions(indices(2), :);

                % Mutation
                mutant = positions(i, :) + (positionRand1 - positionRand2) .* Levy(dim);
                mutant = max(min(mutant, ub), lb); % Ensure within bounds
            else
                % DE_currentToBest
                % Select two random indices
                indices = randperm(searchAgentsNum, 2);
                positionRand1 = positions(indices(1), :);
                positionRand2 = positions(indices(2), :);

                % Mutation
                mutant = positions(i, :) + F * (bestPosition - positions(i, :)) + (positionRand1 - positionRand2) .* Levy(dim);
                mutant = max(min(mutant, ub), lb); % Ensure within bounds

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
            trialFitness = fobj(trial);
            fe = fe + 1;

            if trialFitness < fitness(i)
                positions(i, :) = trial;
                fitness(i) = trialFitness;
            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end
