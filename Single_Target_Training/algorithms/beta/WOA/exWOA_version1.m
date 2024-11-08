% The Whale optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = exWOA_version1(searchAgentsNum, maxFes, lb, ub, dim, fobj)
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

            if abs(A) > 1
                positions(i, :) = positions(i, :) .* Levy(dim) + C * (bestPosition - rand * ones(1, dim));
            elseif abs(A) <= 1

                for j = 1:size(positions, 2)

                    if p < 0.5

                        positionD = abs(C * bestPosition(j) - positions(i, j));
                        positions(i, j) = bestPosition(j) - A * positionD;

                    elseif p >= 0.5
                        distanceToWhale = abs(bestPosition(j) - positions(i, j));
                        positions(i, j) = distanceToWhale * exp(b .* l) .* cos(l .* 2 * pi) + bestPosition(j);

                    end

                end

            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end
