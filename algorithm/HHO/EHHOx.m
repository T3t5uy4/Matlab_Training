% The Enhanced Harris hawks optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = EHHO(searchAgentsNum, maxIters, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    convergenceCurve = [];

    iter = 0;

    while iter < maxIters

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

                end

            end

        end

        % Orthogonal learning
        for i = 1:size(positions, 1)
            bestPosition = orthogonalLearning(positions(i, :), bestPosition, dim);
        end

        % Opposition based learning
        positions = oppositionBasedLearning(positions, dim, fobj);

        for i = 1:size(positions, 1)
            fitness = fobj(positions(i, :));

            if fitness < bestFitness
                bestFitness = fitness;
                bestPosition = positions(i, :);
            end

        end

        iter = iter + 1;
        convergenceCurve(iter) = bestFitness;

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

% Construct an orthogonal table
function L = constructionOfOrthogonalTable(Q, J)
    G = transpose(0:(Q - 1));
    V = zeros(Q ^ J, J);

    for k = 1:J
        R = kron(kron(ones(Q ^ (k - 1), 1), G), ones(Q ^ (J - k), 1));
        V(:, k) = R;
    end

    U = rref(V'); % Simplified row-minimum echelon matrix
    [m, n] = find(U == 1);
    [~, t] = find(U((m + 1):J, n) <= 1);
    B = U(:, unique(n(t)))';
    L = V * B';
    L = mod(L, Q);
    L = L + ones(size(L, 1), size(L, 2));
end

% Orthogonal learning
function bestPosition = orthogonalLearning(position, bestPosition, dim)
    Q = 5;
    J = 2;
    F = (Q ^ J - 1) / (Q - 1);
    L = constructionOfOrthogonalTable(Q, J);
    P = zeros(dim, Q);

    for i = 1:dim
        P(i, 1) = min(position(1, i), rabbitEnergy(1, i));

        for j = 2:(Q - 1)
            P(i, j) = min(position(1, i), rabbitEnergy(1, i)) + (j - 1) * (abs(position(1, i) - bestPosition(1, i)) / (Q - 1));
        end

        P(i, Q) = max(position(1, i), rabbitEnergy(1, i));
    end

end

% Opposition based learning
function positions = oppositionBasedLearning(positions, dim, fobj)
    positionsb = positions;
    k = rand;

    for i = 1:size(positions, 1)
        lb = min(positions(i, :));
        ub = max(positions(i, :));

        for j = 1:dim
            positionsb(i, j) = k * (lb(j) + ub(j)) - positions(i, j);
        end

        if (fobj(positionsb(i, :) < fojb(positions(i, :))))
            positions(i, :) = positionsb(i, :);
        end

    end

end
