% Particle Swarm Optimization - PSO
function [bestFitness, bestPosition, convergenceCurve] = PSO(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    velocities = zeros(searchAgentsNum, dim); % Velocity of each agent
    preBestPosition = positions; % Best position of each agent so far
    fitness = []; % Fitness values for each agent
    convergenceCurve = []; % To record the convergence curve (fitness over iterations)

    w = 0.5; % Inertia weight
    c1 = 1.5; % Cognitive coefficient
    c2 = 1.5; % Social coefficient
    t = 0; % Iteration counter
    fe = 0; % Function evaluations counter

    % Main PSO loop
    while fe < maxFes

        for i = 1:searchAgentsNum
            % Check boundaries and correct positions if out of bounds
            FU = positions(i, :) > ub;
            FL = positions(i, :) < lb;
            positions(i, :) = (positions(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;

            % Evaluate the fitness of the current position
            fitness(i) = fobj(positions(i, :));
            fe = fe + 1;

            % Update the best global position and fitness
            if fitness(i) < bestFitness
                bestFitness = fitness(i);
                bestPosition = positions(i, :);
            end

        end

        % Update particle velocities and positions
        for i = 1:searchAgentsNum
            % Random coefficients
            r1 = rand(dim);
            r2 = rand(dim);

            % Update velocities
            velocities(i, :) = w * velocities(i, :) + c1 * r1 .* (preBestPosition(i, :) - positions(i, :)) + c2 * r2 .* (bestPosition - positions(i, :));

            % Update positions
            positions(i, :) = positions(i, :) + velocities(i, :);

            % Apply bounds
            positions(i, :) = max(min(positions(i, :), ub), lb);

            % Re-evaluate fitness of new positions
            fitness(i) = fobj(positions(i, :));
            fe = fe + 1;

            % Update individual best position
            if fitness(i) < fobj(preBestPosition(i, :))
                preBestPosition(i, :) = positions(i, :);
            end

            fe = fe + 1;

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end
