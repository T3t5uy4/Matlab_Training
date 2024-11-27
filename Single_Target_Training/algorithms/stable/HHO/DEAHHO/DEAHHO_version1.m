% The Harris hawks optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = DEAHHO_version1(searchAgentsNum, maxIters, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    convergenceCurve = [];

    iter = 0;
    t = 0;

    while iter < maxIters

        for i = 1:size(positions, 1)
            % Check boundries
            FU = positions(i, :) > ub;
            FL = positions(i, :) < lb;
            positions(i, :) = (positions(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
            % Fitness of locations
            fitness = fobj(positions(i, :));
            iter = iter + 1;

            if fitness < bestFitness
                bestFitness = fitness;
                bestPosition = positions(i, :);
            end

        end

        % Update the Harris hawk's position
        for i = 1:size(positions, 1)
            % Get the escape energies
            E0 = 2 * rand - 1;
            E = 2 * E0 * (1 - iter / maxIters);

            if abs(E) >= 1 % Exploration stage
                q = rand;
                randHawkIdx = floor(searchAgentsNum * rand + 1);
                randPosition = positions(randHawkIdx, :);

                if q < 0.5
                    % Perch based on other family members
                    positions(i, :) = randPosition - rand * abs(randPosition - 2 * rand * positions(i, :));
                elseif q >= 0.5
                    % Perch on a random tall tree (random site inside group's home range)
                    positions(i, :) = bestPosition - mean(positions) - rand * ((ub - lb) * rand + lb);
                end

            elseif abs(E) < 1 % Development stage
                r = rand;

                if r >= 0.5 && abs(E) < 0.5 % Hard besiege
                    positions(i, :) = bestPosition - E * abs(bestPosition - positions(i, :));
                elseif r >= 0.5 && abs(E) >= 0.5 % Soft besiege
                    jumpStrength = 2 * (1 - rand); % Random jump strength
                    positions(i, :) = bestPosition - positions(i, :) - E * abs(jumpStrength * bestPosition - positions(i, :));
                elseif r < 0.5 && abs(E) >= 0.5 % Soft besiege with motions
                    jumpStrength = 2 * (1 - rand); % Random jump strength
                    position1 = bestPosition - E * abs(jumpStrength * bestPosition - positions(i, :));

                    if fobj(position1) < fobj(positions(i, :))
                        positions(i, :) = position1;
                    else
                        position2 = bestPosition - E * abs(jumpStrength * bestPosition - positions(i, :)) + rand(1, dim) .* Levy(dim);

                        if fobj(position2) < fobj(positions(i, :))
                            positions(i, :) = position2;
                        end

                    end

                    iter = iter + 1;

                elseif r < 0.5 && abs(E) < 0.5 % Hard besiege with motions
                    jumpStrength = 2 * (1 - rand); % Random jump strength
                    position1 = bestPosition - E * abs(jumpStrength * bestPosition - mean(positions));

                    if fobj(position1) < fobj(positions(i, :))
                        positions(i, :) = position1;
                    else
                        position2 = bestPosition - E * abs(jumpStrength * bestPosition - mean(positions)) + rand(1, dim) .* Levy(dim);

                        if fobj(position2) < fobj(positions(i, :))
                            positions(i, :) = position2;
                        end

                    end

                    iter = iter + 1;

                end

            end

        end

        % Sort positions as fobj
        fobjPositions = [];

        for i = 1:size(positions, 1)
            fobjPositions(i) = fobj(positions(i, :));
        end

        [~, sortedIndices] = sort(fobjPositions);
        sortedPositions = positions(sortedIndices, :);

        % Divede the positions
        partSize = floor(searchAgentsNum / 3);

        bestPart = sortedPositions(1:partSize, :);
        worstPart = sortedPositions(2 * partSize + 1:end, :);

        % Update the best part
        [bestPart, iter] = DEA(bestPart, dim, fobj, iter);

        % Update the worst part
        for i = 1:size(worstPart, 1)
            worstPart(i, :) = worstPart(i, :) + 0.05 * (bestPosition - worstPart(i, :));
        end

        % Merge positions
        mergedPositions = [bestPart; positions; worstPart];
        fobjPositions = [];

        for i = 1:size(mergedPositions, 1)
            fobjPositions(i) = fobj(mergedPositions(i, :));
        end

        [~, sortedIndices] = sort(fobjPositions);
        sortedPositions = mergedPositions(sortedIndices, :);
        positions = sortedPositions(1:searchAgentsNum, :);
        fitness = fobj(positions(1, :));

        if fitness < bestFitness
            bestFitness = fitness;
            bestPosition = positions(1, :);
        end

        iter = iter + 1;
        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end

% This function initialize the first population of search agents
function positions = initialization(searchAgentsNum, dim, ub, lb)

    boundaryNum = size(ub, 2); % number of boundaries

    % If the boundaries of all variables are equal and user enter a signle
    % number for both ub and lb
    if boundaryNum == 1
        positions = rand(searchAgentsNum, dim) .* (ub - lb) + lb;
    end

    % If each variable has a different lb and ub
    if boundaryNum > 1

        for i = 1:dim
            ub_i = ub(i);
            lb_i = lb(i);
            positions(:, i) = rand(searchAgentsNum, 1) .* (ub_i - lb_i) + lb_i;
        end

    end

end

% Levi-flight
function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2 ^ ((beta - 1) / 2))) ^ (1 / beta);
    u = randn(1, d) * sigma; v = randn(1, d); step = u ./ abs(v) .^ (1 / beta);
    o = step;
end

% The Differential Evolution Algorithm
function [positions, iter] = DEA(positions, dim, fobj, iter)
    F = 0.5;
    CR = 0.9;
    maxIters = 50;
    newPositions = positions;

    for t = 1:maxIters

        for i = 1:size(newPositions, 1)
            idxs = randperm(size(newPositions, 1), 3);
            r1 = idxs(1);
            r2 = idxs(2);
            r3 = idxs(3);

            mutant = newPositions(r1, :) + F * (newPositions(r2, :) - newPositions(r3, :));

            trial = newPositions(i, :);

            for j = 1:dim

                if rand < CR
                    trial(j) = mutant(j);
                end

            end

            if fobj(trial) < fobj(newPositions(i, :))
                newPositions(i, :) = trial;
                iter = iter + 1;
            end

        end

    end

    tempo = [];

    for i = 1:size(positions, 1)

        if ~isequal(newPositions(i, :), positions(i, :))
            tempo(end + 1, :) = newPositions(i, :);
        end

    end

    positions = tempo;
end
