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
            epsilon = max(epsilon, (fitness(i) - bestFitness) / 2);
        end

        for i = 1:size(positions, 1)
            delta = fitness(i) - bestFitness;
            dp = zeros(1, dim);

            if delta <= 0.2 * (1 - (fe / maxFes)) * epsilon

                for j = 1:dim

                    if bestPosition(1, j) <= positions(i, j)
                        dp(1, j) = 1;
                    else
                        dp(1, j) = -1;
                    end

                end

            elseif delta <= 0.7 * (1 - (fe / maxFes)) * epsilon
                dp = 2 * randi([0, 1], 1, dim) - 1;
            else

                for j = 1:dim

                    if bestPosition(1, j) <= positions(i, j)
                        dp(1, j) = -1;
                    else
                        dp(1, j) = 1;
                    end

                end

            end

            % Adjust the ratio of exploration and exploitation dynamically
            [exploreRatio, exploitRatio] = dynamicPlanning(fe, maxFes);

            % Adjust alpha, beta, and lr dynamically
            [alpha, beta, lr] = dynamicPhaseAdjustment(fe, maxFes);
            r1 = abs(2 * rand - 1) * (1 - fe / maxFes);
            r2 = abs(2 * rand - 1);

            if r1 < exploreRatio
                % Perform global search with random perturbations
                if r2 < 0.5
                    newPosition = lr * positions(i, :) + abs(alpha * ((2 * rand - 1) * (ub - lb) .* rand(1, dim)) + beta * (2 * rand - 1) * bestPosition) .* dp;
                else
                    newPosition = lr * positions(i, :) + abs(alpha * ((2 * rand - 1) * (ub - lb) .* rand(1, dim)) + beta * (2 * rand - 1) * historyBestPositions(i, :)) .* dp;
                end

            elseif r1 < exploreRatio + exploitRatio
                % Perform local search with some global perturbation
                randIdx = getRandIndex(i, searchAgentsNum, 3);
                randPosition1 = positions(randIdx(1), :);
                randPosition2 = positions(randIdx(2), :);

                if r2 <= 1/3
                    newPosition = positions(i, :) + abs(lr * (bestPosition - positions(i, :)) + alpha * ((2 * rand - 1) * (ub - lb) .* rand(1, dim)) + beta * (randPosition1 - randPosition2)) .* dp;
                elseif r2 <= 2/3
                    newPosition = positions(i, :) + abs(lr * (historyBestPositions(i, :) - positions(i, :)) + alpha * ((2 * rand - 1) * (ub - lb) .* rand(1, dim)) + beta * (randPosition1 - randPosition2)) .* dp;
                elseif
                    newPosition = positions(i, :) + abs(lr * (historyBestPositions(randIdx(3), :) - positions(i, :)) + alpha * ((2 * rand - 1) * (ub - lb) .* rand(1, dim)) + beta * (randPosition1 - randPosition2)) .* dp;
                end

            else
                % Perform perturbation with Levy flights if stuck in local optimum
                r3 = rand;
                randIdx = getRandIndex(i, searchAgentsNum, 2);
                randPosition1 = positions(randIdx(1), :);
                randPosition2 = positions(randIdx(2), :);

                if r2 < 0.5
                    newPosition = lr * positions(i, :) + abs(alpha * (randPosition1 - randPosition2) + beta * (2 * rand - 1) * bestPosition + r3 * mean(positions)) .* dp .* Levy(dim);
                else
                    newPosition = lr * positions(i, :) + abs(alpha * (randPosition1 - randPosition2) + beta * (2 * rand - 1) * historyBestPositions(i, :) + r3 * mean(positions)) .* dp .* Levy(dim);
                end

            end

            newFitness = fobj(newPosition);
            fe = fe + 1;

            if newFitness < fitness(i)
                fitness(i) = newFitness;
                positions(i, :) = newPosition;

                if fitness(i) < bestFitness
                    bestFitness = fitness(i);
                    bestPosition = positions(i, :);
                end

                if fitness(i) < historyBestFitness(i)
                    historyBestFitness(i) = fitness(i);
                    historyBestPositions(i, :) = positions(i, :);
                end

            end

            % % The memory mechanism: When the current solution is far from the historical best solution, allow larger jumps in the search space.

            distanceToBest = norm(positions(i, :) - bestPosition);

            if distanceToBest > epsilon % Assume a distance threshold.
                r4 = (exp(1 + (t / maxFes)) + exp(1 - (t / maxFes))) / 2;
            else
                r4 = (exp(1 + (t / maxFes)) - exp(1 - (t / maxFes))) / 2;
            end

            positions(i, :) = positions(i, :) + cos(r4) .* Levy(dim);

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;
    end

end

function [exploreRatio, exploitRatio] = dynamicPlanning(fe, maxFes)
    progress = fe / maxFes;

    if progress <= 1/3
        exploreRatio = 0.7; exploitRatio = 0.2;
    elseif progress <= 2/3
        exploreRatio = 0.3; exploitRatio = 0.4;
    else
        exploreRatio = 0.1; exploitRatio = 0.2;
    end

end

function [lr, alpha, beta] = dynamicPhaseAdjustment(fe, maxFes)
    % Dynamic programming function to adjust parameters based on search phase
    % and historical data (e.g., search progress)

    % Adjust alpha, beta, lr based on current progress and fitness
    progress = fe / maxFes;

    if progress <= 1/3
        lr = 0.5;
        alpha = 0.6;
        beta = 0.8;
    elseif progress <= 2/3
        lr = 0.3;
        alpha = 0.4;
        beta = 0.6;
    else
        lr = 0.1;
        alpha = 0.2;
        beta = 0.4;
    end

end

function [k] = getRandIndex(i, n, count)

    candidates = setdiff(1:n, i);
    k = datasample(candidates, count, 'Replace', false);
end
