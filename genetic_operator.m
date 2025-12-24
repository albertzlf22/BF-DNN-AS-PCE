function f  = genetic_operator(parent_chromosome, grade, popu_num, c_index, m_index, l_limit, u_limit)
% ·±Ö³×Ó´ú
%% function f  = genetic_operator(parent_chromosome, popu_num, c_index, m_index, l_limit, u_limit)
% generate offspring
%
% parent_chromosome - the set of selected chromosomes.
% grade - the non-dominated order and crowded-distance
% popu_num - number of child individuals
% c_index - distribution index for crossover
% m_index - distribution index for mutation
% l_limit - a vector of lower limit for the corresponding decsion variables
% u_limit - a vector of upper limit for the corresponding decsion variables
%
% Perfrom crossover and Mutation operator
% The original NSGA-II uses Simulated Binary Crossover(SBX) and Polynomial mutation.
% Crossover probability pc = 0.9 and mutation probability is pm = 1/n
% where n is the number of decision variables.


%%
[N,V] = size(parent_chromosome);
if mod(popu_num,2)==0
    f=zeros(popu_num,V);
else
    f=zeros(popu_num+1,V);
end
i=0;
while i<popu_num
    child_1=zeros(1,V);
    child_2=zeros(1,V);
    % Select parent
    index=randperm(N,4);
    parent_1=tournament_selection(parent_chromosome(index(1),:),grade(index(1),:),parent_chromosome(index(2),:),grade(index(2),:));
    parent_2=tournament_selection(parent_chromosome(index(3),:),grade(index(3),:),parent_chromosome(index(4),:),grade(index(4),:));
    % Perform corssover
    for j = 1 : V
        % SBX (Simulated Binary Crossover).
        rand_u = rand(1);
        if rand_u <= 0.5
            bq = (2 * rand_u)^(1/(c_index+1));
        else
            bq = (1/(2*(1 - rand_u)))^(1/(c_index+1));
        end
        % Generate the corresponding child element
        child_1(j) = min([max([0.5*(((1 + bq)*parent_1(j)) + (1 - bq)*parent_2(j)),l_limit(j)]),u_limit(j)]);
        child_2(j) = min([max([0.5*(((1 - bq)*parent_1(j)) + (1 + bq)*parent_2(j)),l_limit(j)]),u_limit(j)]);
        %With probability of 1/V doing mutation
        if rand(1)<=1/V
            delta_1=(child_1(j)-l_limit(j))/(u_limit(j)-l_limit(j));
            delta_2=(u_limit(j)-child_1(j))/(u_limit(j)-l_limit(j));
            r=rand(1);
            if r<=0.5
                delta=(2*r+(1-2*r)*(1-delta_1)^(m_index+1))^(1/(m_index+1));
            else
                delta=1-(2*(1-r)+(2*r-1)*(1-delta_2)^(m_index+1))^(1/(m_index+1));
            end
            child_1(j) = min([max([child_1(j)+delta*(u_limit(j)-l_limit(j)),l_limit(j)]),u_limit(j)]);
        end
        if rand(1)<=1/V
            delta_1=(child_2(j)-l_limit(j))/(u_limit(j)-l_limit(j));
            delta_2=(u_limit(j)-child_2(j))/(u_limit(j)-l_limit(j));
            r=rand(1);
            if r<=0.5
                delta=(2*r+(1-2*r)*(1-delta_1)^(m_index+1))^(1/(m_index+1));
            else
                delta=1-(2*(1-r)+(2*r-1)*(1-delta_2)^(m_index+1))^(1/(m_index+1));
            end
            child_2(j) = min([max([child_2(j)+delta*(u_limit(j)-l_limit(j)),l_limit(j)]),u_limit(j)]);
        end
    end
    f(i+1,:)=child_1;
    f(i+2,:)=child_2;
    i=i+2;
end

%% Tournament Selection
    function winner=tournament_selection(player_1,grade_1,player_2,grade_2)
        if grade_1(1) < grade_2(1)
            winner = player_1;
        elseif grade_1(1) > grade_2(1)
            winner = player_2;
        elseif grade_1(2) > grade_2(2)
            winner = player_1;
        elseif grade_1(2) < grade_2(2)
            winner = player_2;
        elseif rand(1)<0.5
            winner = player_1;
        else
            winner = player_2;
        end
    end
end