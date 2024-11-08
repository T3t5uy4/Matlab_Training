% The Harris hawks optimization Algorithm
function [bestFitness, bestPosition, convergenceCurve] = RLHHO(searchAgentsNum, maxFes, lb, ub, dim, fobj)
    % Initialize position vector and fitness for the best
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    % Initialize the positions of search agents
    positions = initialization(searchAgentsNum, dim, ub, lb);
    fitness = [];
    convergenceCurve = [];
    qTable1 = zeros(2, 3);
    rTable1 = [-1 1 1];
    qTable2 = zeros(3, 2);
    rTable2 = [-1 1; 1 1; 1 1];

    t = 0;
    fe = 0;
    F = 0;
    proB = 0.8;
    betaMin = 0.2;
    betaMax = 0.8;
    lbd = 0.1;
    sgm = 0.9;

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
                    fitness1 = fobj(position1);
                    fe = fe + 1;

                    if fitness1 < fitness(i)
                        fitness(i) = fitness1;
                        positions(i, :) = position1;
                    else
                        position2 = bestPosition - E * abs(jumpStrength * bestPosition - positions(i, :)) + rand(1, dim) .* Levy(dim);
                        fitness2 = fobj(position2);
                        fe = fe + 1;

                        if fitness2 < fitness(i)
                            fitness(i) = fitness2;
                            positions(i, :) = position2;
                        end

                    end

                elseif r < 0.5 && abs(E) < 0.5 % Hard besiege with motions
                    jumpStrength = 2 * (1 - rand); % Random jump strength
                    position1 = bestPosition - E * abs(jumpStrength * bestPosition - mean(positions));
                    fitness1 = fobj(position1);
                    fe = fe + 1;

                    if fitness1 < fitness(i)
                        fitness(i) = fitness1;
                        positions(i, :) = position1;
                    else
                        position2 = bestPosition - E * abs(jumpStrength * bestPosition - mean(positions)) + rand(1, dim) .* Levy(dim);
                        fitness2 = fobj(position2);
                        fe = fe + 1;

                        if fitness2 < fitness(i)
                            fitness(i) = fitness2;
                            positions(i, :) = position2;
                        end

                    end

                end

            end

            [k1, k2] = get2Rand(i, searchAgentsNum);

            if rand < proB
                % Update positions with mutation stragey based on improved Q-Learning

                % selecting the valuable mutation operator

                for u = 1:size(qTable1, 1)

                    if qTable1(1, u) >= qTable1(1, 1) && qTable1(1, u) >= qTable1(1, 2) && qTable1(1, u) >= qTable1(1, 3)

                        if u == 1
                            F = normrnd(0.5, 0.3);
                        elseif u == 2
                            F = rand;
                        elseif u == 3
                            F = unifrnd(betaMin, betaMax);
                        end

                        for v = 1:size(qTable2, 2)

                            if qTable2(u, v) >= qTable2(u, 1) && qTable2(u, v) >= qTable2(u, 2)

                                if v == 1
                                    vPosition = bestPosition + F .* (positions(i, :) - positions(k1, :));
                                else
                                    vPosition = bestPosition + F .* (positions(k1, :) - positions(k2, :));
                                end

                                for j = 1:dim
                                    k = sigmf(vPosition(j), [10 0.5]);

                                    if k < rand
                                        k = 0;
                                    else
                                        k = 1;
                                    end

                                    vPosition(j) = k;
                                end

                                vFitness = fobj(vPosition);

                                if vFitness < fitness(i)
                                    fitness(i) = vFitness;
                                    positions(i, :) = vPosition;
                                    rTable1(1, u) = 1;
                                    rTable2(u, v) = 1;
                                else
                                    rTable1(1, u) = -1;
                                    rTable2(u, v) = -1;
                                end

                                tempMax1 = -inf;

                                for c = 1:size(qTable1, 2)

                                    if qTable1(1, c) > tempMax1
                                        tempMax1 = qTable1(1, c);
                                    end

                                end

                                tempMax2 = -inf;

                                for c = 1:size(qTable2, 2)

                                    if qTable2(1, c) > tempMax1
                                        tempMax1 = qTable2(1, c);
                                    end

                                end

                                qTable1(1, u) = qTable1(1, u) + lbd * (rTable1(1, u) + sgm * tempMax1 - qTable1(1, u));
                                qTable2(u, v) = qTable2(u, v) + lbd * (rTable2(u, v) + sgm * tempMax2 - qTable2(u, v));
                                fe = fe + 1;

                                break;
                            end

                        end

                        break;
                    end

                end

            else
                % Update positions with randomly selected mutation stragegy
                u = randi([1, 3]);
                v = randi([1, 2]);

                if u == 1
                    F = normrnd(0.5, 0.3);
                elseif u == 2
                    F = rand;
                elseif u == 3
                    F = unifrnd(betaMin, betaMax);
                end

                if v == 1
                    vPosition = bestPosition + F .* (positions(i, :) - positions(k1, :));
                else
                    vPosition = bestPosition + F .* (positions(k1, :) - positions(k2, :));
                end

                for j = 1:dim
                    k = sigmf(vPosition(j), [10 0.5]);

                    if k < rand
                        k = 0;
                    else
                        k = 1;
                    end

                    vPosition(j) = k;
                end

                vFitness = fobj(vPosition);

                if vFitness < fitness(i)
                    fitness(i) = vFitness;
                    positions(i, :) = vPosition;
                    rTable1(1, u) = 1;
                    rTable2(u, v) = 1;
                else
                    rTable1(1, u) = -1;
                    rTable2(u, v) = -1;
                end

                tempMax1 = -inf;

                for c = 1:size(qTable1, 2)

                    if qTable1(1, c) > tempMax1
                        tempMax1 = qTable1(1, c);
                    end

                end

                tempMax2 = -inf;

                for c = 1:size(qTable2, 2)

                    if qTable2(1, c) > tempMax1
                        tempMax1 = qTable2(1, c);
                    end

                end

                qTable1(1, u) = qTable1(1, u) + lbd * (rTable1(1, u) + sgm * tempMax1 - qTable1(1, u));
                qTable2(u, v) = qTable2(u, v) + lbd * (rTable2(u, v) + sgm * tempMax2 - qTable2(u, v));
                fe = fe + 1;

            end

        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;

    end

end

function [k1, k2] = get2Rand(i, n)
    k1 = randi([1, n]);

    while k1 == i
        k1 = randi([1, n]);
    end

    k2 = randi([1, n]);

    while k2 == i || k2 == k1
        k2 = randi([1, n]);
    end

end
