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
    epsilon = (ub - lb) / exp(dim / 2);

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

        % Adjust the ratio of exploration and exploitation dynamically
        [exploreRatio, exploitRatio] = dynamicPlanning(fe, maxFes);

        % Adjust alpha, beta, and lr dynamically
        [lr, alpha, beta] = dynamicPhaseAdjustment(fe, maxFes);

        for i = 1:size(positions, 1)
            r1 = rand * (fe / maxFes);
            r2 = rand;

            if r1 < exploreRatio
                % Perform global search with random perturbations
                if r2 < 0.5
                    newPosition = positions(i, :) + lr * (bestPosition - positions(i, :)) + alpha * (lb + (ub - lb) .* rand(1, dim)) + beta .* Levy(dim);
                else
                    newPosition = positions(i, :) + lr * (historyBestPositions(i, :) - positions(i, :)) + alpha * (lb + (ub - lb) .* rand(1, dim)) + beta .* Levy(dim);
                end

            elseif r1 < exploreRatio + exploitRatio
                % Perform local search with some global perturbation
                randIdx = getRandIndex(i, searchAgentsNum, 2);
                randPosition1 = positions(randIdx(1), :);
                randPosition2 = positions(randIdx(2), :);

                if r2 < 0.5
                    newPosition = positions(i, :) + lr * (bestPosition - positions(i, :)) + alpha * (lb + (ub - lb) .* rand(1, dim)) + beta * (randPosition1 - randPosition2);
                else
                    newPosition = positions(i, :) + lr * (historyBestPositions(i, :) - positions(i, :)) + alpha * (lb + (ub - lb) .* rand(1, dim)) + beta * (randPosition1 - randPosition2);
                end

            else
                % Perform perturbation with Levy flights if stuck in local optimum
                r3 = 2 * rand - 1;

                if r2 < 0.5
                    newPosition = bestPosition + lr * (historyBestPositions(i, :) - positions(i, :)) + alpha * (lb + (ub - lb) .* rand(1, dim)) + beta * r3 * positions(i, :);
                else
                    newPosition = historyBestPositions(i, :) + lr * (bestPosition - positions(i, :)) + alpha * (lb + (ub - lb) .* rand(1, dim)) + beta * r3 * positions(i, :);
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

            % The memory mechanism: When the current solution is far from the historical best solution, allow larger jumps in the search space.
            distanceToBest = norm(positions(i, :) - bestPosition);

            if distanceToBest > epsilon % Assume a distance threshold.
                r4 = 2 * rand - 1;
                positions(i, :) = positions(i, :) + r4 .* Levy(dim);
            else
                positions(i, :) = positions(i, :) + sqrt(2) .* randn(1, dim);
            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;
    end

end

function [exploreRatio, exploitRatio] = dynamicPlanning(fe, maxFes)
    progress = fe / maxFes;

    if progress <= 1/3
        exploreRatio = 0.8; exploitRatio = 0.15;
    elseif progress <= 2/3
        exploreRatio = 0.5; exploitRatio = 0.4;
    else
        exploreRatio = 0.3; exploitRatio = 0.6;
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
        beta = 0.9;
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
