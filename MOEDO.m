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
clc;
clear;
close all;
% Problem Configuration
ObjectiveFunction = @ZDT1; % Objective function handle
dim = 2; % Number of dimensions
ub = ones(1, dim); % Upper bounds (assuming 1 for all dimensions)
lb = zeros(1, dim); % Lower bounds (assuming 0 for all dimensions)
obj_no = 2; % Number of objectives
% Algorithm Parameters
max_iter = 100; % Maximum number of iterations
ArchiveMaxSize = 100; % Maximum size of the archive
Archive_X = zeros(ArchiveMaxSize, dim); % Initialize archive solutions
Archive_F = ones(ArchiveMaxSize, obj_no) * inf; % Initialize archive fitnesses
Archive_member_no = 0; % Number of members in the archive
Xwinners = initialization(ArchiveMaxSize, dim, ub, lb); % Population initialization
Fitness = zeros(ArchiveMaxSize, obj_no); % obj_no is the number of objectives, which is 2
for i = 1:ArchiveMaxSize
    Fitness(i, :) = ObjectiveFunction(Xwinners(i,:));
end
Memoryless=Xwinners;
% Main loop for MOEDO algorithm
for iter = 1:max_iter
    V=zeros(ArchiveMaxSize,dim);
    [Fitness,sorted_indices]=sort(Fitness);
    temp_Xwinners=Xwinners;
    Xwinners=temp_Xwinners(sorted_indices,:);
    d=(1-iter/max_iter);  
    f= 2*rand-1; 
    a=f^10;      
    b=f^5;       
    c=d*f;       
    sum=(Xwinners(1,:)+Xwinners(2,:)+Xwinners(3,:));
    X_guide=(sum/3);  
    for  i=1:ArchiveMaxSize
        alpha=rand;
        if alpha<0.5
            if Memoryless(i,:)==Xwinners(i,:)
                Mu=(X_guide+Memoryless(i,:))/2.0;   
                ExP_rate=1./Mu;    
                variance=1./ExP_rate.^2;   
                V(i,:)=a.*(Memoryless(i,:)-variance)+b.*X_guide;   
            else
                Mu=(X_guide+Memoryless(i,:))/2.0;    
                ExP_rate=1./Mu;   
                variance=1./ExP_rate.^2;   
                phi=rand;
                V(i,:)=b.*(Memoryless(i,:)-variance)+log (phi).*Xwinners(i,:);
            end
        else
            M=mean(Xwinners);
            s=randperm(ArchiveMaxSize);
            D1=M-Xwinners(s(1),:); 
            D2=M-Xwinners(s(2),:); 
            Z1=M-D1+D2; 
            Z2=M-D2+D1; 
            V(i,:)=(Xwinners(i,:)+(c.*Z1+(1-c).*Z2)-M); 
        end
        F_UB=V(i,:)>ub;
        F_LB=V(i,:)<lb;
        V(i,:)=(V(i,:).*(~(F_UB+F_LB)))+ub.*F_UB+lb.*F_LB;
    end
    for i=1:ArchiveMaxSize
        for j=1:dim
            Memoryless(i,j)=V(i,j);
        end
    end
    V_Fitness = zeros(ArchiveMaxSize, obj_no); % obj_no is the number of objectives
for i = 1:ArchiveMaxSize
    V_Fitness(i, :) = ObjectiveFunction(V(i,:));
        Xwinners(i, :) = V(i, :);
        Fitness(i, :) = V_Fitness(i, :);
end
    % Non-dominated Sorting and Crowding Distance Calculation
    Combined_X = [Xwinners; Archive_X(1:Archive_member_no, :)];
    Combined_F = [V_Fitness; Archive_F(1:Archive_member_no, :)];
    [fronts, ~] = NonDominatedSorting(Combined_F);
    crowdingDistances = CrowdingDistance(Combined_F, fronts);
    % Update Archive using NSGA-II strategies
    [Archive_X, Archive_F, Archive_member_no] = UpdateArchiveUsingNSGAII(fronts, crowdingDistances, Combined_X, Combined_F, ArchiveMaxSize);
    % Display iteration information
    disp(['At iteration ', num2str(iter), ', MOEDO has ', num2str(Archive_member_no), ' non-dominated solutions in the archive']);
end
% Plotting the results
figure;
Draw_ZDT1(); % Function to draw the true Pareto Front (assuming it is defined)
hold on;
plot(Archive_F(:, 1), Archive_F(:, 2), 'ro', 'MarkerSize', 8, 'markerfacecolor', 'k');
legend('True PF', 'Obtained PF');
title('MOEDO');
set(gcf, 'pos', [403 466 230 200]); % Setting the figure position and size
