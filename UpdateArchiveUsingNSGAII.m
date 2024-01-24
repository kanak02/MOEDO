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


function [Archive_X, Archive_F, Archive_member_no] = UpdateArchiveUsingNSGAII(fronts, distances, Combined_X, Combined_F, ArchiveMaxSize)
    Archive_X = [];
    Archive_F = [];
    Archive_member_no = 0;
    for f = 1:length(fronts)
        front = fronts{f};
        [~, I] = sort(-distances(front));
        front = front(I);
        if length(front) + Archive_member_no <= ArchiveMaxSize
            Archive_X = [Archive_X; Combined_X(front, :)];
            Archive_F = [Archive_F; Combined_F(front, :)];
            Archive_member_no = Archive_member_no + length(front);
        else
            remaining = ArchiveMaxSize - Archive_member_no;
            Archive_X = [Archive_X; Combined_X(front(1:remaining), :)];
            Archive_F = [Archive_F; Combined_F(front(1:remaining), :)];
            Archive_member_no = ArchiveMaxSize;
            break;
        end
    end
end
