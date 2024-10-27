% Particle Swarm Optimization - PSO
function [bestFitness, bestPosition, convergenceCurve] = PSO(searchAgentNum, maxFes, lb, ub, dim, fobj)

    % Define the PSO's paramters
    ub = ub .* ones(1, dim);
    lb = lb .* ones(1, dim);

    wMax = 0.9;
    wMin = 0.2;
    vMax = (ub - lb) * 0.2;
    vMin = -vMax;
    c1 = 2;
    c2 = 2;

    % Initializations
    iter = maxFes;
    pBestScore = zeros(searchAgentNum);
    pBest = zeros(searchAgentNum, dim);
    bestPosition = zeros(1, dim);
    convergenceCurve = zeros(1, iter);
    velocity = zeros(searchAgentNum, dim);
    positions = zeros(searchAgentNum, dim);

    %Initialization
    for i = 1:size(positions, 1)

        for j = 1:size(positions, 2)
            positions(i, j) = (ub(j) - lb(j)) * rand() + lb(j);
            velocity(i, j) = rand();
        end

    end

    for i = 1:searchAgentNum
        pBestScore(i) = inf;
    end

    % Initialize gBestScore for a minimization problem
    bestFitness = inf;

    for t = 1:iter

        for i = 1:size(positions, 1)
            % Return back the particles that go beyond the boundaries of the search
            % space
            Flag4ub = positions(i, :) > ub;
            Flag4lb = positions(i, :) < lb;
            positions(i, :) = (positions(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            %Calculate objective function for each particle
            fitness = fobj(positions(i, :));

            if (pBestScore(i) > fitness)
                pBestScore(i) = fitness;
                pBest(i, :) = positions(i, :);
            end

            if (bestFitness > fitness)
                bestFitness = fitness;
                bestPosition = positions(i, :);
            end

        end

        %Update the W of PSO
        w = wMax - t * ((wMax - wMin) / iter);
        %Update the Velocity and Position of particles
        for i = 1:size(positions, 1)

            for j = 1:size(positions, 2)
                velocity(i, j) = w * velocity(i, j) + c1 * rand() * (pBest(i, j) - positions(i, j)) + c2 * rand() * (bestPosition(j) - positions(i, j));

                if velocity(i, j) > vMax(j)
                    velocity(i, j) = vMax(j);
                end

                if velocity(i, j) <- vMax(j)
                    velocity(i, j) = -vMax(j);
                end

                positions(i, j) = positions(i, j) + velocity(i, j);
            end

        end

        convergenceCurve(t) = bestFitness;
    end

end
