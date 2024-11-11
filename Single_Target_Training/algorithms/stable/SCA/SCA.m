% The SCA optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = SCA(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    fitness = [];
    convergenceCurve = [];

    t = 0;
    fe = 0;

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

        a = 2;
        r1 = a - t * (a / maxFes);

        % 更新位置
        for i = 1:size(positions, 1)

            for j = 1:size(positions, 2)
                r2 = 2 * pi * rand();
                r3 = 2 * rand();
                r4 = rand();

                if r4 < 0.5
                    positions(i, j) = positions(i, j) + r1 * sin(r2) * abs(r3 * bestPosition(j) - positions(i, j));
                else
                    positions(i, j) = positions(i, j) + r1 * cos(r2) * abs(r3 * bestPosition(j) - positions(i, j));
                end

            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end
