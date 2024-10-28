% The Whale optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = MWOA(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    fitness = [];
    convergenceCurve = [];

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

        a = 2 * (1 - fe / maxFes);
        a2 = -1 - (fe / maxFes);
        % Update the Whale's position
        for i = 1:size(positions, 1)
            r1 = rand;
            r2 = rand;

            A = 2 * a * r1 - a;
            C = 2 * r2;

            b = 1;
            l = (a2 - 1) * rand + 1;
            p = rand;

            for j = 1:size(positions, 2)

                if p < 0.5

                    if abs(A) >= 1
                        randWhaleIdx = floor(searchAgentsNum * rand + 1);
                        randPosition = positions(randWhaleIdx, :);
                        randPositionD = abs(C * randPosition(j) - positions(i, j));
                        positions(i, j) = randPosition(j) - A * randPositionD;
                    elseif abs(A) < 1
                        positionD = abs(C * bestPosition(j) - positions(i, j));
                        positions(i, j) = bestPosition(j) - A * positionD;
                    end

                elseif p >= 0.5
                    distanceToWhale = abs(bestPosition(j) - positions(i, j));
                    positions(i, j) = distanceToWhale * exp(b .* l) .* cos(l .* 2 * pi) + bestPosition(j);
                end

            end

        end

        % Levy-flight
        for i = 1:size(positions, 1)
            positionV = positions(i, :) + 0.75 * Levy(dim) .* (positions(i, :) - bestPosition);
            fitnessV = fobj(positionV);
            fe = fe + 1;

            if fitnessV < fobj(positions(i, :))
                positions(i, :) = positionV;
            end

            if fitnessV < bestFitness
                bestFitness = fitnessV;
                bestPosition = positionV;
            end

        end

        % pattern search
        if mod(fe, round(maxFes * 0.1)) == 0
            options = optimoptions('patternsearch', 'Display', 'off');
            [positionPS, fitnessPS] = patternsearch(fobj, bestPosition, [], [], [], [], [], [], [], options);

            if fitnessPS < bestFitness
                bestFitness = fitnessPS;
                bestPosition = positionPS;
            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end
