function [bestFitness, bestPosition, Convergence_curve] = RUN(searchAgentNum, maxFes, lb, ub, dim, fobj)

    Cost = zeros(searchAgentNum, 1); % Record the Fitness of all Solutions
    X = initialization(searchAgentNum, dim, ub, lb); % Initialize the set of random solutions
    Xnew2 = zeros(1, dim);

    Convergence_curve = zeros(1, maxFes);

    for i = 1:searchAgentNum
        Cost(i) = fobj(X(i, :)); % Calculate the Value of Objective Function
    end

    [bestFitness, ind] = min(Cost); % Determine the Best Solution
    bestPosition = X(ind, :);

    Convergence_curve(1) = bestFitness;

    %% Main Loop of RUN
    it = 1; %Number of iterations

    while it < maxFes
        it = it + 1;
        f = 20 .* exp(- (12 .* (it / maxFes))); % (Eq.17.6)
        Xavg = mean(X); % Determine the Average of Solutions
        SF = 2 .* (0.5 - rand(1, searchAgentNum)) .* f; % Determine the Adaptive Factor (Eq.17.5)

        for i = 1:searchAgentNum
            [~, ind_l] = min(Cost);
            lBest = X(ind_l, :);

            [A, B, C] = RndX(searchAgentNum, i); % Determine Three Random Indices of Solutions
            [~, ind1] = min(Cost([A B C]));

            % Determine Delta X (Eqs. 11.1 to 11.3)
            gama = rand .* (X(i, :) - rand(1, dim) .* (ub - lb)) .* exp(-4 * it / maxFes);
            Stp = rand(1, dim) .* ((bestPosition - rand .* Xavg) + gama);
            DelX = 2 * rand(1, dim) .* (abs(Stp));

            % Determine Xb and Xw for using in Runge Kutta method
            if Cost(i) < Cost(ind1)
                Xb = X(i, :);
                Xw = X(ind1, :);
            else
                Xb = X(ind1, :);
                Xw = X(i, :);
            end

            SM = RungeKutta(Xb, Xw, DelX); % Search Mechanism (SM) of RUN based on Runge Kutta Method

            L = rand(1, dim) < 0.5;
            Xc = L .* X(i, :) + (1 - L) .* X(A, :); % (Eq. 17.3)
            Xm = L .* bestPosition + (1 - L) .* lBest; % (Eq. 17.4)

            vec = [1, -1];
            flag = floor(2 * rand(1, dim) + 1);
            r = vec(flag); % An Interger number

            g = 2 * rand;
            mu = 0.5 + .1 * randn(1, dim);

            % Determine New Solution Based on Runge Kutta Method (Eq.18)
            if rand < 0.5
                Xnew = (Xc + r .* SF(i) .* g .* Xc) + SF(i) .* (SM) + mu .* (Xm - Xc);
            else
                Xnew = (Xm + r .* SF(i) .* g .* Xm) + SF(i) .* (SM) + mu .* (X(A, :) - X(B, :));
            end

            % Check if solutions go outside the search space and bring them back
            FU = Xnew > ub; FL = Xnew < lb; Xnew = (Xnew .* (~(FU + FL))) + ub .* FU + lb .* FL;
            CostNew = fobj(Xnew);

            if CostNew < Cost(i)
                X(i, :) = Xnew;
                Cost(i) = CostNew;
            end

            %% Enhanced solution quality (ESQ)  (Eq. 19)
            if rand < 0.5
                EXP = exp(-5 * rand * it / maxFes);
                r = floor(Unifrnd(-1, 2, 1, 1));

                u = 2 * rand(1, dim);
                w = Unifrnd(0, 2, 1, dim) .* EXP; %(Eq.19-1)

                [A, B, C] = RndX(searchAgentNum, i);
                Xavg = (X(A, :) + X(B, :) + X(C, :)) / 3; %(Eq.19-2)

                beta = rand(1, dim);
                Xnew1 = beta .* (bestPosition) + (1 - beta) .* (Xavg); %(Eq.19-3)

                for j = 1:dim

                    if w(j) < 1
                        Xnew2(j) = Xnew1(j) + r * w(j) * abs((Xnew1(j) - Xavg(j)) + randn);
                    else
                        Xnew2(j) = (Xnew1(j) - Xavg(j)) + r * w(j) * abs((u(j) .* Xnew1(j) - Xavg(j)) + randn);
                    end

                end

                FU = Xnew2 > ub; FL = Xnew2 < lb; Xnew2 = (Xnew2 .* (~(FU + FL))) + ub .* FU + lb .* FL;
                CostNew = fobj(Xnew2);

                if CostNew < Cost(i)
                    X(i, :) = Xnew2;
                    Cost(i) = CostNew;
                else

                    if rand < w(randi(dim))
                        SM = RungeKutta(X(i, :), Xnew2, DelX);
                        Xnew = (Xnew2 - rand .* Xnew2) + SF(i) * (SM + (2 * rand(1, dim) .* bestPosition - Xnew2)); % (Eq. 20)

                        FU = Xnew > ub; FL = Xnew < lb; Xnew = (Xnew .* (~(FU + FL))) + ub .* FU + lb .* FL;
                        CostNew = fobj(Xnew);

                        if CostNew < Cost(i)
                            X(i, :) = Xnew;
                            Cost(i) = CostNew;
                        end

                    end

                end

            end

            % End of ESQ
            %% Determine the Best Solution
            if Cost(i) < bestFitness
                bestPosition = X(i, :);
                bestFitness = Cost(i);
            end

        end

        % Save Best Solution at each iteration
        Convergence_curve(it) = bestFitness;
        disp(['it : ' num2str(it) ', Best Cost = ' num2str(Convergence_curve(it))]);

    end

end

% A funtion to determine a random number
%with uniform distribution (unifrnd function in Matlab)
function z = Unifrnd(a, b, c, dim)
    a2 = a / 2;
    b2 = b / 2;
    mu = a2 + b2;
    sig = b2 - a2;
    z = mu + sig .* (2 * rand(c, dim) - 1);
end

% A function to determine thress random indices of solutions
function [A, B, C] = RndX(nP, i)
    Qi = randperm(nP); Qi(Qi == i) = [];
    A = Qi(1); B = Qi(2); C = Qi(3);
end

function SM = RungeKutta(XB, XW, DelX)

    dim = size(XB, 2);
    C = randi([1 2]) * (1 - rand);
    r1 = rand(1, dim);
    r2 = rand(1, dim);

    K1 = 0.5 * (rand * XW - C .* XB);
    K2 = 0.5 * (rand * (XW + r2 .* K1 .* DelX / 2) - (C * XB + r1 .* K1 .* DelX / 2));
    K3 = 0.5 * (rand * (XW + r2 .* K2 .* DelX / 2) - (C * XB + r1 .* K2 .* DelX / 2));
    K4 = 0.5 * (rand * (XW + r2 .* K3 .* DelX) - (C * XB + r1 .* K3 .* DelX));

    XRK = (K1 + 2 .* K2 + 2 .* K3 + K4);
    SM = 1/6 * XRK;
end
