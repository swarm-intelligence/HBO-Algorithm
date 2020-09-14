clc;
clear;
close all;

%Number of Initial Population & Dimensions
pop_size = 100;
dim = 2;

%Iteration Condition
max_iter = 10;

%Domain of Benchmarks
from = -32;
to = -1*from;

%Results for n Times Execution
num_of_result = 5;
%Columns of total Result : dim,(gbest_fitness),(time)
total_result = zeros(num_of_result,dim+2);

for q=1:num_of_result
    tic;
    
    %Nfe Condition
    max_nfe = 20000;
    
    %Initilize Best Fitness and Position by a Large Value
    best = zeros(1,dim);
    best(:,:) = to;
    best_fit = F2(best); %nfe++
    nfe = 1;
    
    %Determine the Size of m Patches from n
    m_size = round(0.75*pop_size);
    %Determine the Size of Elite Patches
    e_size = round(0.75*m_size);
    
    %Main Loop
    for a=1:max_iter
        %Initial Population
        n = unifrnd(from, to, [pop_size dim]);
        e = zeros(e_size,dim);
        s = size(n,1);
        
        %Calculate Fitness of n and Sorting
        [m, nfe, F_result_l1] = HBO_Sort(n,s,m_size,nfe,dim,max_nfe);
        
        %Select the Elite Patches
        e(:,:) = m(1:e_size,:);
        %The Remaining Population from m (m-e)
        m_e = m(e_size+1:end,:);
        
        %Determine the Number of NEP and NSP Bees
        n_nep = ceil((size(n,1)/(size(e,1)+size(m,1))) + size(e,1));
        n_nep = n_nep * size(e,1);
        n_nsp = ceil(size(m_e,1)/2);
        n_nsp = n_nsp * size(m_e,1);
        
        %Calculate Neighborhood Radius for each Iteration
        ngh = 1 + a/(max_iter+a);
        %alpha0 = 0.5 + (1-0.5)*rand;
        %ngh = alpha0 + a/max_iter;
        
        nep = zeros(n_nep,dim);
        nsp = zeros(n_nsp,dim);
        
        %Scatter the NEP Bees Around the e Patches
        for i=1:size(e,1)
            for j=1:dim
                %Creating the NEP Bees with NGH
                nep(i,j) = e(i,j)-ngh + ((e(i,j)+ngh) - (e(i,j)-ngh))*rand;
                
                s = size(nep,1);
                %Calculate the Fintness of NEP Bees
                [nep, nfe, F_result_nep] = HBO_Sort(nep,s,n_nep,nfe,dim,max_nfe);
                
                if(nfe >= max_nfe)
                    break;
                end
                
                %Select the Best NEP
                best_nep = F_result_nep(1);
                
                %Compare the Best NEP with the Best Previous Answer
                if(best_nep < best_fit)
                    best_fit = best_nep;
                    best(1,:) = nep(1,:);
                end
            end
        end
        
        %Scatter the NSP Bees Around the m-e Patches
        for i=1:size(m_e,1)
            for j=1:dim
                %Creating the NSP Bees with NGH
                nsp(i,j) = m_e(i,j)-ngh + ((m_e(i,j)+ngh) - (m_e(i,j)-ngh))*rand;
                
                s = size(nsp,1);
                %Calculate the Fintness of NSP Bees
                [nsp, nfe, F_result_nsp] = HBO_Sort(nsp,s,n_nsp,nfe,dim,max_nfe);
                
                if(nfe >= max_nfe)
                    break;
                end
                
                %Select the Best NSP
                best_nsp = F_result_nsp(1);
                
                %Compare the Best NEP with the Best Previous Answer
                if(best_nsp < best_fit)
                    best_fit = best_nsp;
                    best(1,:) = nsp(1,:);
                end
            end
        end
        
        %Global Search (Random Search)
        for c=1:pop_size
            %Determine the Size of Random Answers
            random_pop = ceil(0.1*pop_size);
            %Create the Random Answers
            global_random = unifrnd(from, to, [random_pop dim]);
            
            s = size(global_random,1);
            %Calculate the Fitness of Random Answers
            [global_random, nfe, F_result_g] = HBO_Sort(global_random,s,random_pop,nfe,dim,max_nfe);
            
            if(nfe >= max_nfe)
                break;
            end
            
            %Select the Best Random Answer
            best_global = F_result_g(1);
            
            %Compare the Best Random Answer with the Best Previous Answer
            if(best_global < best_fit)
                best_fit = F_result_g;
                best(1,:) = global_random(1,:);
            end
            
        end
    end
    
    total_result(q,1) = toc;
    total_result(q,2) = best_fit;
    total_result(q,3:end) = best;
end

min_fitness = min(total_result(:,2));
max_fitness = max(total_result(:,2));
mean_fitness = mean(total_result(:,2));
std_fitness = std(total_result(:,2));
mean_time = mean(total_result(:,1));

disp(strcat('Popsize:', num2str(pop_size), ', Dim:', num2str(dim)));
disp(strcat('mean fitness: ', num2str(mean_fitness)));
disp(strcat('max fitness: ', num2str(max_fitness)));
disp(strcat('min fitness: ', num2str(min_fitness)));
disp(strcat('std fitness: ', num2str(std_fitness)));
disp(strcat('mean time: ', num2str(mean_time)));