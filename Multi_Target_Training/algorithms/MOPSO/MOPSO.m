function [bestFitness, bestPosition, paretoFront] = MOPSO(searchAgentsNum, maxFes, lb, ub, dim, fobj)

    positions = lb + (ub - lb) * rand(searchAgentsNum, dim);
    velocities = zeros(searchAgentsNum, dim);
    paretoFront = [];
    t = 0;
    fe = 0;

    while fe < maxFes

    end

end
