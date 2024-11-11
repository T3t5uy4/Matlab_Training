% The Harris hawks optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = CXSHHO(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);
    secondFitness = inf;
    secondPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    fitness = [];
    convergenceCurve = [];

    t = 0;
    fe = 0;
    alpha = 0.95;

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
                secondFitness = bestFitness;
                bestFitness = fitness(i);
                secondPosition = bestPosition;
                bestPosition = positions(i, :);
            elseif fitness(i) < secondFitness
                secondFitness = fitness(i);
                secondPosition = positions(i, :);
            end

        end

        for i = 1:size(positions, 1)
            % Communication and collaboration strategy
            r6 = rand;
            beta = r6 / alpha;
            X(i, :) = X(i, :) + (1 - r6) * (bestPosition - secondPosition) / 2;
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

            % Directional crossover strategy
            [k1, k2] = getRand(i, searchAgentsNum);
            parent1 = positions(k1, :);
            parent2 = positions(k2, :);

            for j = 1:size(positions, 2)

                if parent1(j) ~= parent2(j)

                    if fobj(parent1(j)) <= fobj(parent2(j))
                        bestParent = parent1(j);
                    else
                        bestParent = parent2(j)
                    end

                    meanParent = (parent1(j) + parent2(j)) / 2;

                    if bestParent >= meanParent
                        val = 1 - 0.5 ^ (exp(abs(parent1(j) - parent2(j)) / (bestPosition(j) - secondPosition(j))));
                        c1 = val * (parent1(j) - parent2(j)) + alpha ^ (1 - r6) * exp(1 - beta) * (1 - val) * abs(parent1(j) - parent2(j));
                        c2 = (1 - val) * (parent1(j) - parent2(j));

                    else
                    end

                else
                end

            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end

function [k1, k2] = get2Rand(i, n)
    k1 = randi([1, n]);

    while k1 == i
        k1 = randi([1, n]);
    end

    k2 = randi([1, n]);

    while k2 == i || k2 == k1
        k2 = randi([1, n]);
    end

end
