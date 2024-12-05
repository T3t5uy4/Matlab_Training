% Tortoise and Hare optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = THO(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    fitness = inf * ones(1, dim);
    convergenceCurve = [];

    t = 0;
    fe = 0;
    r = (ub - lb) / 2 * rand;
    omega = 0.1;
    k = 0.1;
    theta = rand(1, dim) * 2 * pi;
    alpha = 0.6;
    beta = 0.2;
    gamma = 0.1;

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

        % Reselecting the hare and the tortoise
        [fitness, idx] = sort(fitness);

        if fe / maxFes > 1/2
            % Exchanging the order of precedence between the tortoise and the hare
            mid = floor(size(positions, 1) / 2);
            idx = [idx(mid + 1:end), idx(1:mid)];
        end

        positions = positions(idx, :);

        theta = theta + omega;
        r = r + k * (fe / maxFes);

        % Update the location of the hare
        for i = 1:size(positions, 1) / 2

            if fe / maxFes <= 1/5
                % Hare Dumps Tortoise Stage
                positions(i, :) = positions(i, :) + alpha * rand(1, dim) .* Levy(dim);
            elseif fe / maxFes <= 4/5
                % Lazy rabbit stage
                positions(i, :) = positions(i, :) + r .* (cos(theta) - 0.5);
            else
                % Hare Chasing Tortoise Stage
                distance = norm(bestPosition - positions(i, :));
                positions(i, :) = positions(i, :) + beta * (1 - fe / maxFes) * (bestPosition - positions(i, :)) / (distance + 0.01);
            end

        end

        % Update the location of the tortoise
        for i = size(positions, 1) / 2 + 1:size(positions, 1)

            if fe / maxFes <= 1/2
                % Tortoise chasing hare stage
                q = rand;

                if q <= 0.5
                    positions(i, :) = positions(i, :) + gamma * (bestPosition - positions(i, :));
                else
                    hareIdx = randi([1, size(positions, 1) / 2]);
                    positions(i, :) = positions(i, :) + gamma * (positions(hareIdx, :) - positions(i, :));
                end

            else
                % Tortoise overtakes the hare stage
                positions(i, :) = positions(i, :) + (2 * rand(1, dim) - 1) .* rand(1, dim);
            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end
