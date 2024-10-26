% The Particle Swarm Optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = PSO(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    convergenceCurve = [];

    fe = 0;

    while fe < maxFes

        for i = 1:size(positions, 1)
            % Check boundries
            FU = positions(i, :) > ub;
            FL = positions(i, :) < lb;
            positions(i, :) = (positions(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
            % Fitness of locations
            fitness = fobj(positions(i, :));

            if fitness < bestFitness
                bestFitness = fitness;
                bestPosition = positions(i, :);
            end

        end

        % Update

        fe = fe + 1;
        convergenceCurve(iter) = bestFitness;

    end

end
