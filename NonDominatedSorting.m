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
function [fronts, maxFront] = NonDominatedSorting(F)
    % Initialize
    [S, n, frontNumbers] = deal(cell(size(F, 1), 1));
    [rank, distances] = deal(zeros(size(F, 1), 1));
    front = 1;
    maxFront = 0;

    % Calculate domination
    for i = 1:size(F, 1)
        S{i} = [];
        n{i} = 0;
        for j = 1:size(F, 1)
            if dominates(F(i, :), F(j, :))
                S{i} = [S{i}, j];
            elseif dominates(F(j, :), F(i, :))
                n{i} = n{i} + 1;
            end
        end
        if n{i} == 0
            rank(i) = 1;
            if isempty(frontNumbers{front})
                frontNumbers{front} = i;
            else
                frontNumbers{front} = [frontNumbers{front}, i];
            end
        end
    end

    % Assign fronts
    while ~isempty(frontNumbers{front})
        Q = [];
        for i = frontNumbers{front}
            for j = S{i}
                n{j} = n{j} - 1;
                if n{j} == 0
                    rank(j) = front + 1;
                    Q = [Q, j];
                end
            end
        end
        front = front + 1;
        frontNumbers{front} = Q;
    end
    maxFront = front - 1;

    % Organize fronts
    fronts = cell(maxFront, 1);
    for i = 1:maxFront
        fronts{i} = frontNumbers{i};
    end
end
