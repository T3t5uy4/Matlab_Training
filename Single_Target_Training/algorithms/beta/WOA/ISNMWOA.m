% The Whale Optimization Algorithm
function [Leader_score, Leader_pos, Convergence_curve] = ISNMWOA(SearchAgents_no, MaxFEs, lb, ub, dim, fobj)
    Leader_pos = zeros(1, dim);
    Leader_score = inf;
    J = 1/3;
    alpha = 0.5;
    Positions = initialization(SearchAgents_no, dim, ub, lb);
    lb = ones(1, dim) .* lb;
    ub = ones(1, dim) .* ub;
    Convergence_curve = [];
    FEs = 0;
    t = 1;

    while FEs < MaxFEs

        for i = 1:size(Positions, 1)
            Flag4ub = Positions(i, :) > ub;
            Flag4lb = Positions(i, :) < lb;
            Positions(i, :) = (Positions(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            fitness = fobj(Positions(i, :));
            FEs = FEs + 1;

            if fitness < Leader_score
                Leader_score = fitness;
                Leader_pos = Positions(i, :);
            end

        end

        for i = 1:round(SearchAgents_no)
            p = rand;
            indexJ = randi(SearchAgents_no);

            if p < J
                Positions(i, :) = Positions(i, :) +0.01 .* (Positions(indexJ, :) .* Levy(dim) - Positions(i, :));
            elseif p < (1 - J)
                Positions(i, :) = Positions(i, :) + alpha .* (Positions(indexJ, :) - Positions(i, :) .* Levy(dim));
            else
                Positions(i, :) = Positions(i, :) + alpha .* (Positions(indexJ, :) - Positions(i, :)) .* Levy(dim);
            end

        end

        a = 2 - FEs * ((2) / MaxFEs);
        a2 = -1 + FEs * ((-1) / MaxFEs);

        for i = 1:size(Positions, 1)
            r1 = rand();
            r2 = rand();
            A = 2 * a * r1 - a;
            C = 2 * r2;
            b = 1;
            l = (a2 - 1) * rand + 1;
            p = rand();

            for j = 1:size(Positions, 2)

                if p < 0.5

                    if abs(A) >= 1
                        rand_leader_index = floor(SearchAgents_no * rand() + 1);
                        X_rand = Positions(rand_leader_index, :);
                        D_X_rand = abs(C * X_rand(j) - Positions(i, j));
                        Positions(i, j) = X_rand(j) - A * D_X_rand;
                    elseif abs(A) < 1
                        D_Leader = abs(C * Leader_pos(j) - Positions(i, j));
                        Positions(i, j) = Leader_pos(j) - A * D_Leader;
                    end

                elseif p >= 0.5
                    distance2Leader = abs(Leader_pos(j) - Positions(i, j));
                    Positions(i, j) = distance2Leader * exp(b .* l) .* cos(l .* 2 * pi) + Leader_pos(j);
                end

            end

        end

        options = optimset('MaxFunEvals', floor(MaxFEs * 0.1));
        [x, fval, ~, output] = fminsearchbnd(fobj, Leader_pos, lb, ub, options);

        if fval < Leader_score
            Leader_score = fval;
            Leader_pos = x;
        end

        FEs = FEs + output.funcCount;
        Convergence_curve(t) = Leader_score;
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
