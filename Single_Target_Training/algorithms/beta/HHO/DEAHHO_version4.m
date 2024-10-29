% The Harris hawks optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = DEAHHO_version4(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    convergenceCurve = [];

    fe = 0;
    t = 0;
    keepRate = 0.0;

    while fe < maxFes

        for i = 1:size(positions, 1)
            % Check boundries
            FU = positions(i, :) > ub;
            FL = positions(i, :) < lb;
            positions(i, :) = (positions(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
            % Fitness of locations
            fitness = fobj(positions(i, :));
            fe = fe + 1;

            if fitness < bestFitness
                bestFitness = fitness;
                bestPosition = positions(i, :);
            end

        end

        % Update the Harris hawk's position
        for i = 1:size(positions, 1)
            % Get the escape energies
            E0 = 2 * rand - 1;
            E = 2 * E0 * (1 - fe / maxFes);

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

                    fe = fe + 1;

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

                    fe = fe + 1;

                end

            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end

% The Differential Evolution Algorithm
function [positions, fe] = DEA(positions, dim, fobj, fe)
    F = 0.5;
    CR = 0.9;
    maxFes = 50;
    newPositions = positions;

    for t = 1:maxFes

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
                fe = fe + 1;
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
