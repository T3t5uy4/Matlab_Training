function [bestFitness, bestPosition, convergenceCurve] = no_name(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);
    secondFitness = inf;
    secondPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    convergenceCurve = []
    fitness = zeros(searchAgentsNum, 1);
    t = 0;
    fe = searchAgentsNum;
    alpha = 0.6;
    beta = 0.2;
    lr = 0.1;
    preCount = 0;
    learningTable = [2 1];

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
                secondPosition = bestPosition;
                bestFitness = fitness(i);
                bestPosition = positions(i, :);
            elseif fitness(i) < secondFitness
                secondFitness = fitness(i);
                secondPosition = positions(i, :);
            end

        end

        for i = 1:size(positions, 1)
            p = fe / maxFes;
            r1 = rand;
            r2 = rand;
            count = 0;

            if r1 > p
                positions(i, :) = alpha * positions(i, :) + learningTable(1) * beta * (rand(1, dim) .* (ub - lb) + lb .* ones(1, dim)) + learningTable(2) * lr * (bestPosition + secondPosition) / 2;
            else if r1 <= p
                indices = getRandIndex(i, searchAgentsNum, 1);
                positions(i, :) = alpha * positions(i, :) + learningTable(1) * beta * (rand(1, dim) .* (ub - lb) + lb .* ones(1, dim)) + learningTable(2) * lr * positions(indices(1), :);
            end

            fitness(i) = fobj(positions(i, :));
            fe = fe + 1;

            if fitness(i) < bestFitness
                secondFitness = bestFitness;
                secondPosition = bestPosition;
                bestFitness = fitness(i);
                bestPosition = positions(i, :);
                count = count + 2;
            elseif fitness(i) < secondFitness
                secondFitness = fitness(i);
                secondPosition = positions(i, :);
            end

        end

        if count > preCount
            learningTable(1) = 1;
        elseif count < preCount

            if learningTable(1) == 2
                learningTable(2) = -1 * learningTable(2);
            end

            learningTable(1) = 2;
        end

        for i = 1:size(positions, 1)
            r2 = rand;

            if r2 >= count / (searchAgentsNum * 2)
                positions(i, :) = positions(i, :) .* Levy(dim);
            else
                r3 = rand;
                indices = getRandIndex(i, searchAgentsNum, 4);
                positionRand1 = positions(indices(1), :);
                positionRand2 = positions(indices(2), :);
                positionRand3 = positions(indices(3), :);
                positionRand4 = positions(indices(4), :);
                positions(i, :) = positions(i, :) + r3 * (positionRand1 - positionRand2) + (1 - r3) * (positionRand3 - positionRand4);
            end

        end

        preCount = count;
        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end

function [k] = getRandIndex(idx, n, count)

    for i = 1:count
        k(i) = rand([1, n]);

        while k(i) == i
            k(i) = rand([1, n]);
        end

    end

end
