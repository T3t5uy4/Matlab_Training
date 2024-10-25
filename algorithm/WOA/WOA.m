% The Whale optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = WOA(searchAgentsNum, maxIters, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    convergenceCurve = [];

    iter = 0;

    while iter < maxIters

        for i = 1:size(positions, 1)
            % Check boundries
            FU = positions(i, :) > ub;
            FL = positions(i, :) < lb;
            positions(i, :) = (positions(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
            % Fitness of locations
            fitness = fobj(positions(i, :));

            if fitness < bestFitness
                bestFitness = fitness;
                bestPosition = positions(i, :);
            end

        end

        a = 2 * (1 - iter / maxIters);
        a2 = -1 - (iter / maxIters);
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

        iter = iter + 1;
        convergenceCurve(iter) = bestFitness;

    end

end

% This function initialize the first population of search agents
function positions = initialization(searchAgentsNum, dim, ub, lb)

    boundaryNum = size(ub, 2); % number of boundaries

    % If the boundaries of all variables are equal and user enter a signle
    % number for both ub and lb
    if boundaryNum == 1
        positions = rand(searchAgentsNum, dim) .* (ub - lb) + lb;
    end

    % If each variable has a different lb and ub
    if boundaryNum > 1

        for i = 1:dim
            ub_i = ub(i);
            lb_i = lb(i);
            positions(:, i) = rand(searchAgentsNum, 1) .* (ub_i - lb_i) + lb_i;
        end

    end

end
