% The Harris hawks optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = RLHHO(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    fitness = [];
    convergenceCurve = [];
    qTable1 = zeros(3, 2);
    rTable1 = [-1 1 1];
    QTable2 = zeros(2, 3);
    rTable2 = [-1 1; 1 1; 1 1];

    t = 0;
    fe = 0;
    proB = 0.8;
    betaMin = 0.2;
    betaMax = 0.8;
    lbd = 0.1;
    sgm = 0.9;

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

        end

        if rand < proB
            % Update positions with mutation stragey based on improved Q-Learning

            % selecting the valuable mutation operator
            if qTable11(1, 1) >= qTable1(1, 2) && qTable1(1, 1) >= qTable1(1, 3)
                F = normrnd(0.5, 0.3);

                if qTable2(1, 1) >= qTable2(1, 2)

                end

            elseif qTable1(1, 2) >= qTable1(1, 1) && qTable1(1, 2) >= QTable2(1, 3)
                F = rand;
            else
                F = unifrnd(betaMin, betaMax, dim);
            end

        else
            % Update positions with randomly selected mutation stragegy
        end

        % Update Q-Table

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end
