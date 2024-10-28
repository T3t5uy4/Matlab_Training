% The Whale Optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = ISNMWOA(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    bestPosition = zeros(1, dim);
    bestFitness = inf;
    J = 1/3;
    alpha = 0.5;
    positions = initialization(searchAgentsNum, dim, ub, lb);
    lb = ones(1, dim) .* lb;
    ub = ones(1, dim) .* ub;
    convergenceCurve = [];
    fitness = [];
    fe = 0;
    t = 1;

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

        for i = 1:round(searchAgentsNum)
            p = rand;
            indexJ = randi(searchAgentsNum);

            if p < J
                positions(i, :) = positions(i, :) +0.01 .* (positions(indexJ, :) .* Levy(dim) - positions(i, :));
            elseif p < (1 - J)
                positions(i, :) = positions(i, :) + alpha .* (positions(indexJ, :) - positions(i, :) .* Levy(dim));
            else
                positions(i, :) = positions(i, :) + alpha .* (positions(indexJ, :) - positions(i, :)) .* Levy(dim);
            end

        end

        a = 2 - fe * ((2) / maxFes);
        a2 = -1 + fe * ((-1) / maxFes);

        for i = 1:size(positions, 1)
            r1 = rand();
            r2 = rand();
            A = 2 * a * r1 - a;
            C = 2 * r2;
            b = 1;
            l = (a2 - 1) * rand + 1;
            p = rand();

            for j = 1:size(positions, 2)

                if p < 0.5

                    if abs(A) >= 1
                        rand_leader_index = floor(searchAgentsNum * rand() + 1);
                        X_rand = positions(rand_leader_index, :);
                        D_X_rand = abs(C * X_rand(j) - positions(i, j));
                        positions(i, j) = X_rand(j) - A * D_X_rand;
                    elseif abs(A) < 1
                        D_Leader = abs(C * bestPosition(j) - positions(i, j));
                        positions(i, j) = bestPosition(j) - A * D_Leader;
                    end

                elseif p >= 0.5
                    distance2Leader = abs(bestPosition(j) - positions(i, j));
                    positions(i, j) = distance2Leader * exp(b .* l) .* cos(l .* 2 * pi) + bestPosition(j);
                end

            end

        end

        options = optimset('MaxFunEvals', floor(maxFes * 0.1));
        [x, fval, ~, output] = fminsearchbnd(fobj, bestPosition, lb, ub, options);

        if fval < bestFitness
            bestFitness = fval;
            bestPosition = x;
        end

        fe = fe + output.funcCount;
        convergenceCurve(t) = bestFitness;
        t = t + 1;
    end

end

% Levy-flight
function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2 ^ ((beta - 1) / 2))) ^ (1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v) .^ (1 / beta);
    o = 0.01 * step;
end
