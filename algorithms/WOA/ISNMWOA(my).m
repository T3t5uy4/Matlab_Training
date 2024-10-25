% The Whale optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = ISNMWOA(searchAgentsNum, maxIters, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    convergenceCurve = [];
    lb = ones(1, dim) .* lb;
    ub = ones(1, dim) .* ub;

    iter = 0;
    J = 1/3;

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

        for i = 1:size(positions, 1)
            p = rand;
            randWhaleIdx = floor(searchAgentsNum * rand + 1);

            if p < J
                positions(i, :) = positions(i, :) + 0.01 .* (positions(randWhaleIdx, :) .* Levy(dim) - positions(i, :));
            elseif p < (1 - J)
                positions(i, :) = positions(i, :) + 0.5 .* (positions(randWhaleIdx, :) - positions(i, :) .* Levy(dim));
            else
                positions(i, :) = positions(i, :) + 0.5 .* (positions(randWhaleIdx, :) - positions(i, :)) .* Levy(dim);
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

        % options = optimset('', 'final');
        [bestPosition, bestFitness] = fminsearchbnd(fobj, bestPosition, lb, ub);

        iter = iter + 1;
        convergenceCurve(iter) = bestFitness;

    end

end
