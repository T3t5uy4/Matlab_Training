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

        theta = theta + omega;
        r = r + k * (fe / maxFes);
        alpha = 0.6 * (1 - (fe / maxFes));
        beta = 0.2 * (1 - (fe / maxFes));
        gamma = 0.1 * (1 - (fe / maxFes));
        harSize = floor(size(positions, 1) * max(0.1, (1 - (fe / maxFes))));
        torSize = searchAgentsNum - harSize;
        % Reselecting the hare and the tortoise
        [fitness, idx] = sort(fitness);

        if fe / maxFes <= 1/2
            % Exchanging the order of precedence between the tortoise and the hare
            idx = [idx(1:torSize), idx(torSize + 1:end)];
        end

        positions = positions(idx, :);

        % Update the location of the hare
        for i = 1:harSize

            if fe / maxFes <= 1/5
                % Hare Dumps Tortoise Stage
                positions(i, :) = positions(i, :) + alpha * rand(1, dim) .* Levy(dim);
            elseif fe / maxFes <= 4/5
                % Lazy rabbit stage
                q = rand;

                if q <= 0.5
                    positions(i, :) = positions(i, :) + r .* (cos(theta) - 0.5);
                else
                    positions(i, :) = positions(i, :) + r .* (sin(theta) - 0.5);
                end

            else
                % Hare Chasing Tortoise Stage
                distance = norm(bestPosition - positions(i, :));
                positions(i, :) = positions(i, :) + beta * (1 - fe / maxFes) * (bestPosition - positions(i, :)) / (distance + 0.01);
            end

        end

        % Update the location of the tortoise
        for i = harSize + 1:size(positions, 1)

            if fe / maxFes <= 1/2
                % Tortoise chasing hare stage
                q = rand;

                if q <= 0.5
                    positions(i, :) = positions(i, :) + gamma * (bestPosition - positions(i, :));
                else
                    hareIdx = getRandIndex(i, harSize, 1);
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

function [k] = getRandIndex(i, n, count)
    k = zeros(1, count);

    for i = 1:count
        k(i) = randi([1, n]);

        while k(i) == i
            k(i) = randi([1, n]);
        end

    end

end
