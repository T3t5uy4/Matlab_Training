% Particle Swarm Optimization - PSO
function [bestFitness, bestPosition, convergenceCurve] = PSO(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    convergenceCurve = [];
    velocities = zeros(searchAgentsNum, dim);
    preBestPosition = positions;
    fitness = zeros(searchAgentsNum, 1);

    w = 0.5; % Inertia weight
    c1 = 1.5; % Cognitive coefficient
    c2 = 1.5; % Social coefficient
    t = 0;
    fe = searchAgentsNum;

    while fe < maxFes

        for i = 1:size(positions, 1)
            % Check boundries
            FU = positions(i, :) > ub;
            FL = positions(i, :) < lb;
            positions(i, :) = (positions(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
            % Fitness of locations
            newFitness = fobj(positions(i, :));
            fe = fe + 1;

            if newFitness < bestFitness
                bestFitness = newFitness;
                bestPosition = positions(i, :);
            end

        end

        for i = 1:size(positions, 1)
            % Update velocities
            r1 = rand(dim);
            r2 = rand(dim);
            velocities(i, :) = w * velocities(i, :) + c1 * r1 .* (preBestPosition(i, :) - positions(i, :)) + c2 * r2 .* (bestPosition - positions(i, :));

            % Update positions
            positions(i, :) = positions(i, :) + velocities(i, :);

            % Apply bounds
            positions(i, :) = max(min(positions(i, :), ub), lb);

            % Fitness of locations
            fitness = fobj(positions(i, :));
            fe = fe + 1;

            if fitness < fobj(preBestPosition(i, :))

            end

            if fitness < bestFitness
                bestFitness = fitness;
                bestPosition = positions(i, :);
            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end
