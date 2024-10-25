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
