%____________________________________________________________________________________
% Multi-objective exponential distribution optimizer (MOEDO) source codes version 1.0 (NDS+CD)
%  Author and programmer: Pradeep Jangir
%  Authors:- Kanak Kalita, Janjhyam Venkata Naga Ramesh, Lenka Cepova, Sundaram B. Pandya, Pradeep Jangir & Laith Abualigah
%         e-Mail: pkjmtech@gmail.com
%   Multi-objective exponential distribution optimizer (MOEDO): a novel math-inspired multi-objective algorithm for global optimization and real-world engineering design problems
%   Kanak Kalita, Janjhyam Venkata Naga Ramesh, Lenka Cepova, Sundaram B. Pandya, Pradeep Jangir & Laith Abualigah
%   Scientific Reports volume 14, Article number: 1816 (2024)
%   DOI:https://doi.org/10.1038/s41598-024-52083-7%
%____________________________________________________________________________________

function Positions = initialization(SearchAgents_no, dim, ub, lb)
    Boundary_no = size(ub, 2); % Number of boundaries

    % If the boundaries of all variables are equal and user enters a single
    % number for both ub and lb
    if Boundary_no == 1
        ub_new = ones(1, dim) * ub;
        lb_new = ones(1, dim) * lb;
    else
        ub_new = ub;
        lb_new = lb;   
    end

    % If each variable has a different lb and ub
    Positions = zeros(SearchAgents_no, dim);
    for i = 1:dim
        ub_i = ub_new(i);
        lb_i = lb_new(i);
        Positions(:, i) = rand(SearchAgents_no, 1) .* (ub_i - lb_i) + lb_i;
    end

    % Removed the transpose to ensure correct dimensions
end
