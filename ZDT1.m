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



% Modify this file with respect to your objective function
function o = ZDT1(x)

o = [0, 0];

dim = length(x);
g = 1 + 9*sum(x(2:dim))/(dim-1);

o(1) = x(1);
o(2) = g*(1-sqrt(x(1)/g));



