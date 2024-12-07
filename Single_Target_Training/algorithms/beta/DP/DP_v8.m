function [bestFitness, bestPosition, convergenceCurve] = DP(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    bestFitness = inf;
    bestPosition = zeros(1, dim);
    positions = initialization(searchAgentsNum, dim, ub, lb);
    convergenceCurve = [];
    fitness = zeros(searchAgentsNum, 1);
    historyBestFitness = inf * ones(searchAgentsNum, 1);
    historyBestPositions = positions;
    t = 0;
    fe = 0;

    while fe < maxFes

        for i = 1:size(positions, 1)
            FU = positions(i, :) > ub;
            FL = positions(i, :) < lb;
            positions(i, :) = (positions(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
            fitness(i) = fobj(positions(i, :));
            fe = fe + 1;

            if fitness(i) < bestFitness
                bestFitness = fitness(i);
                bestPosition = positions(i, :);
            end

            if fitness(i) < historyBestFitness(i)
                historyBestFitness(i) = fitness(i);
                historyBestPositions(i, :) = positions(i, :);
            end

        end

        epsilon = -inf;

        for i = 1:size(positions, 1)
            epsilon = max(epsilon, (fitness(i) - bestFitness) * (1 - fe / maxFes));
        end

        for i = 1:size(positions, 1)
            delta = fitness(i) - bestFitness;
            dp = zeros(1, dim);

            if delta <= 0.5 * epsilon

                for j = 1:dim

                    if bestPosition(1, j) <= positions(i, j)
                        dp(1, j) = 1;
                    else
                        dp(1, j) = -1;
                    end

                end

            else

                for j = 1:dim

                    if bestPosition(1, j) <= positions(i, j)
                        dp(1, j) = -1;
                    else
                        dp(1, j) = 1;
                    end

                end

            end

            [alpha, beta, gamma] = dynamicPharameterAdjustment(fe, maxFes);
            r = rand;

            if r <= 1/3
                positions(i, :) = alpha * positions(i, :) + beta * bestPosition .* dp;
            elseif r <= 2/3
                positions(i, :) = alpha * positions(i, :) + beta * historyBestPositions(i, :) .* dp;
            else
                idx = getRandIndex(i, searchAgentsNum, 1);
                positions(i, :) = alpha * positions(i, :) + beta * positions(idx, :) .* dp;
            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;
    end

end

function [alpha, beta, gamma] = dynamicPharameterAdjustment(fe, maxFes)
    % Dynamic programming function to adjust parameters based on search phase
    % and historical data (e.g., search progress)

    % Adjust alpha, beta, lr based on current progress and fitness
    progress = fe / maxFes;

    if progress <= 1/3
        alpha = 0.5;
        beta = 0.6;
        gamma = 0.1;
    elseif progress <= 2/3
        alpha = 0.4;
        beta = 0.5;
        gamma = 0.05;
    else
        alpha = 0.6;
        beta = 0.7;
        gamma = 0.01;
    end

end

function [k] = getRandIndex(i, n, count)

    candidates = setdiff(1:n, i);
    k = datasample(candidates, count, 'Replace', false);
end
