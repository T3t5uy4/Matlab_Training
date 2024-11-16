function [bestFitness, bestPosition, convergenceCurve] = no_name_with_ [alpha, beta, lr](searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);
    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    convergenceCurve = [];
    fitness = zeros(searchAgentsNum, 1);
    t = 0;
    fe = 0;

    while fe < maxFes

        for i = 1:size(positions, 1)
            % Check boundaries
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

        % Dynamic parameter adjustment based on iteration progress
        meanPosition = mean(positions);
        [p1, p2, p3] = adjustR(fe / maxFes);

        % Dynamic Programming to optimize parameters and position update
        for i = 1:size(positions, 1)
            r1 = rand;

            % Calculate DP values for better exploration/exploitation balance
            [alpha, beta, lr] = dynamicProgramming(positions(i, :), bestPosition, meanPosition, fobj, fe, maxFes);

            % DP-based position update
            if r1 <= p1
                positions(i, :) = alpha * positions(i, :) + beta * (lb + (ub - lb) .* rand(1, numel(bestPosition))) + lr * (bestPosition - positions(i, :));
            elseif r1 <= p1 + p2
                positions(i, :) = alpha * positions(i, :) + beta * (lb + (ub - lb) .* rand(1, numel(bestPosition))) + lr * ((bestPosition + meanPosition) / 2 - positions(i, :));
            elseif r1 <= p1 + p2 + p3
                r2 = rand;
                positions(i, :) = alpha * positions(i, :) + beta * (lb + (ub - lb) .* rand(1, numel(bestPosition))) - lr * (bestPosition - positions(i, :)) + r2 * Levy(dim);
            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;
    end

end

% Dynamic programming function to adjust alpha, beta, lr dynamically
function [alpha, beta, lr] = dynamicProgramming(currentPosition, bestPosition, meanPosition, fobj, currentFe, maxFes)
    % Example DP function to decide alpha, beta, lr based on the current state
    % and the iteration progress.

    % Define the state variables based on current search phase
    progress = currentFe / maxFes;

    % Dynamic adjustment of parameters based on the current iteration progress
    if progress <= 1/3
        alpha = 0.8;
        beta = 0.5;
        lr = 0.05;
    elseif progress <= 2/3
        alpha = 0.6;
        beta = 0.3;
        lr = 0.05;
    else
        alpha = 0.4;
        beta = 0.1;
        lr = 0.1;
    end

    % Adjustments based on proximity to global best and mean position
    % Example of a dynamic adjustment: the closer the position to the global best, the smaller the lr (learning rate).
    distanceToBest = norm(currentPosition - bestPosition);
    lr = lr * (1 + exp(-distanceToBest)); % Dynamic scaling of lr

    % Return the DP values
end

function [p1, p2, p3] = adjustR(t)
    % Define the nodes
    t0 = 0; % Initial iteration (start)
    t1 = 0.5; % Midway point (middle)
    t2 = 1; % Final iteration (end)

    % Define the values for each node
    p1_0 = 0.6; p2_0 = 0.2; p3_0 = 0.2;
    p1_1 = 0.2; p2_1 = 0.6; p3_1 = 0.2;
    p1_2 = 0.2; p2_2 = 0.2; p3_2 = 0.6;

    % Calculate p1, p2, p3 using Lagrange interpolation
    p1 = p1_0 * ((t - t1) * (t - t2)) / ((t0 - t1) * (t0 - t2)) + ...
        p1_1 * ((t - t0) * (t - t2)) / ((t1 - t0) * (t1 - t2)) + ...
        p1_2 * ((t - t0) * (t - t1)) / ((t2 - t0) * (t2 - t1));

    p2 = p2_0 * ((t - t1) * (t - t2)) / ((t0 - t1) * (t0 - t2)) + ...
        p2_1 * ((t - t0) * (t - t2)) / ((t1 - t0) * (t1 - t2)) + ...
        p2_2 * ((t - t0) * (t - t1)) / ((t2 - t0) * (t2 - t1));

    p3 = p3_0 * ((t - t1) * (t - t2)) / ((t0 - t1) * (t0 - t2)) + ...
        p3_1 * ((t - t0) * (t - t2)) / ((t1 - t0) * (t1 - t2)) + ...
        p3_2 * ((t - t0) * (t - t1)) / ((t2 - t0) * (t2 - t1));
end
