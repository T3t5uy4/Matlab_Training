% The Rotating Gyro Optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = RGO(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = uniformInitialization(searchAgentsNum, dim, ub, lb);
    fitness = inf * ones(1, dim);
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

        if fe * 2 <= maxFes

            for i = 1:size(positions, 1)
                R = max(positions(i, :)) - min(positions(i, :)) + rand;
                center = mean(positions(i, :));

                for j = 1:dim
                    r = rand;

                    if r <= 0.5
                        positions(i, j) = center +cos(2 * pi * (fe / maxFes)) * R * positions(i, j);
                    else
                        positions(i, j) = center +sin(2 * pi * (fe / maxFes)) * R * positions(i, j);
                    end

                end

            end

        else

            for i = 1:size(positions, 1)

                idx = getRandIndex(i, searchAgentsNum, 3);
                alpha = rand;
                positions(i, :) = cos(2 * pi * positions(i, :)) .* positions(i, :) + positions(idx(1), :) .* Levy(dim) + alpha * (positions(idx(2), :) - positions(idx(3), :));

            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end

function X = uniformInitialization(searchAgentsNum, dim, ub, lb)
    % uniformInitialization - Initializes search agents uniformly in the search space.
    %
    % Syntax: X = uniformInitialization(searchAgentsNum, dim, ub, lb)
    %
    % Inputs:
    %   searchAgentsNum - Number of search agents.
    %   dim             - Dimensionality of the search space.
    %   ub              - Upper bound (scalar or vector of size [1, dim]).
    %   lb              - Lower bound (scalar or vector of size [1, dim]).
    %
    % Outputs:
    %   X - Initialized positions of search agents (matrix of size [searchAgentsNum, dim]).

    % Ensure ub and lb are vectors of length dim
    if isscalar(ub)
        ub = repmat(ub, 1, dim);
    end

    if isscalar(lb)
        lb = repmat(lb, 1, dim);
    end

    % Initialize the search agent positions matrix
    X = zeros(searchAgentsNum, dim);

    % Calculate the uniform interval for each dimension
    stepSizes = (ub - lb) / (searchAgentsNum - 1);

    % Generate evenly spaced values in each dimension
    for i = 1:dim
        X(:, i) = linspace(lb(i), ub(i), searchAgentsNum)';
    end

    % If the number of particles is more than the number of dimensions, choose the right combination of particle positions
    if searchAgentsNum > dim
        X = repmat(X, ceil(searchAgentsNum / dim), 1); % repeat
        X = X(1:searchAgentsNum, :); % keep search agents number
    end

end

function [k] = getRandIndex(i, n, count)
    k = zeros(1, count);

    for i = 1:count
        k(i) = randi([1, n]);

        while k(i) == i
            k(i) = randi([1, n]);
        end

    end

end

function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2 ^ ((beta - 1) / 2))) ^ (1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v) .^ (1 / beta);
    o = 0.01 * step;
end
