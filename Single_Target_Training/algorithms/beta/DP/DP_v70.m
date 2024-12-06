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
    worstFitness = -inf;

    while fe < maxFes

        for i = 1:size(positions, 1)
            FU = positions(i, :) > ub;
            FL = positions(i, :) < lb;
            positions(i, :) = (positions(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
            fitness(i) = fobj(positions(i, :));
            worstFitness = max(worstFitness, fitness(i));
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

            r = rand;

            if fe * 2 <= maxFes

                if r <= 0.7
                    positions(i, :) = positions(i, :) + (2 * rand) * abs(bestPosition - positions(i, :)) .* dp;
                else
                    idx = getRandIndex(i, searchAgentsNum, 1);
                    positions(i, :) = positions(i, :) + (2 * rand) * (positions(idx, :) - positions(i, :)) .* Levy(dim);
                end

            else

                if r <= 0.3
                    positions(i, :) = positions(i, :) + (1 - fe / maxFes) * abs(bestPosition - positions(i, :)) .* dp;
                else
                    idx = getRandIndex(i, searchAgentsNum, 1);
                    positions(i, :) = positions(i, :) + (2 * rand) * (positions(idx, :) - positions(i, :)) .* Levy(dim);
                end

            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;
    end

end

function [k] = getRandIndex(i, n, count)

    candidates = setdiff(1:n, i);
    k = datasample(candidates, count, 'Replace', false);
end
