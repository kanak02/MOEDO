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

function distances = CrowdingDistance(F, fronts)
    % Initialize the distances array
    distances = zeros(size(F, 1), 1);

    % Calculate crowding distance for each front
    for f = 1:length(fronts)
        front = fronts{f};
        frontSize = length(front);

        % Set the boundary points' distances to infinity
        if frontSize > 2
            distances(front) = 0;
            for m = 1:size(F, 2) % Iterate over each objective
                [sortedValues, sortOrder] = sort(F(front, m));
                sortedFront = front(sortOrder);

                % Distance for boundary points
                distances(sortedFront(1)) = inf;
                distances(sortedFront(end)) = inf;

                % Distance for intermediate points
                for i = 2:(frontSize - 1)
                    distances(sortedFront(i)) = distances(sortedFront(i)) + ...
                        (F(sortedFront(i + 1), m) - F(sortedFront(i - 1), m)) / ...
                        (max(F(:, m)) - min(F(:, m)));
                end
            end
        else
            distances(front) = inf;
        end
    end
end
