function f = non_domination_sort(chromosome, popu_num, M, V)
% 非支配排序
%% function f = non_domination_sort_mod(x, M, V)
% Non-Dominated sort
%
% M - inputing the number of objectives
% V - numebr of decision variables
%
% The fast sort algorithm [1]
% 1.Intialize
%   For each individual p in main population P
%     individual.n - the number of individuals that dominate p.被支配数
%     individual.p - all the individuals that is being dominated by p.所支配人口
%
% 2.Judge(minimize)
%   For each individual q in P
%      if p dominated q (minimize)
%         add q to the set i.p i.e. i.p = i.p + {q}.
%      else if q dominates p
%         increment the domination counter i.e. i.n = i.n + 1.
%
% 3.Sort
%   i = 1
%   For each individual p in P
%     if i.n = 0
%         adding p to front one i.e F1 = F1 + {p}.
%         set rank of p to one i.e p.rank = 1.
%
%   While the ith front is nonempty i.e. Fi != []
%     Q = [].
%     For each individual p in front Fi
%         For each individual q in Pp
%            i.n = i.n - 1.
%            if i.n = 0.
%               update the set Q with individual q i.e. Q = Q + {q}.
%               set rank of p to one i.e p.rank = i + 1.
%      i = i + 1.
%      Now the set Q is the next front and hence Fi = Q.
%

%%
N = size(chromosome,1);
% % used only to manipulate easily in MATLAB.
% F(front).f = [];
% individual = [];
dominate=zeros(N,N); %支配标记矩阵，(i,j)==1表示i支配j

for i = 1 : N
    for j = i+1 : N
        % 得到Pareto支配
%         dom_less = 0;
%         dom_equal = 0;
%         dom_more = 0;
%         for m = 1 : M
%             if (chromosome(i,V + m) < chromosome(j,V + m))
%                 dom_less = dom_less + 1;
%             elseif (chromosome(i,V + m) == chromosome(j,V + m))
%                 dom_equal = dom_equal + 1;
%             else
%                 dom_more = dom_more + 1;
%             end
%         end
        dom_less = sum(chromosome(i,V+1:V+M) < chromosome(j,V+1:V+M));
        dom_equal = sum(chromosome(i,V+1:V+M) == chromosome(j,V+1:V+M));
        dom_more = M - dom_less - dom_equal;
        if dom_less == 0 && dom_equal ~= M
            dominate(j,i)=1;
        elseif dom_more == 0 && dom_equal ~= M
            dominate(i,j)=1;
        end
    end
end

chosen_num=0;
front=1;
% order_number = M + V + 1;
% distance_number = M + V + 2;
f=zeros(popu_num,M + V + 2);
% chromosome(:,order_number)=zeros(N,1);
% non-dominated sort
while chosen_num<popu_num
    front_index=find(sum(dominate)==0);    %当前front索引
    front_num=length(front_index);         %当前front个体数量
    front_i=chromosome(front_index,:);%当前front染色体与目标值
    dominate(front_index,:)=zeros(front_num,N);%去除当前front对其他个体的支配标记
    dominate(front_index,front_index)=dominate(front_index,front_index)+diag(inf.*ones(1,front_num));%将当前front在支配标记矩阵的对角元素(i,i)置为inf，表示已被选择
%% 拥挤度计算
    distance=zeros(front_num,1);
    for i = 1 : M
        [c_sorted_on_objective, index_of_objectives] = sort(front_i(:,V + i));
        obj_min = front_i(index_of_objectives(1),V + i);
        obj_max = front_i(index_of_objectives(front_num),V + i);
        if obj_min==obj_max
            distance=inf.*ones(front_num,1);
            break
        end
        distance(index_of_objectives(1)) = inf;        %最小值对应拥挤度设为inf
        distance(index_of_objectives(front_num)) = inf;%最大值对应拥挤度设为inf
        for j = 2 : front_num-1
            next_obj  = c_sorted_on_objective(j + 1);
            previous_obj  = c_sorted_on_objective(j - 1);
            distance(index_of_objectives(j))=distance(index_of_objectives(j))+(next_obj - previous_obj)/(obj_max - obj_min);
        end
    end
%%
%     distance = rand(front_num,1);   %替代拥挤度，随机排序
    [~,index_of_distance] = sort(distance,'descend');
    chosen_num=chosen_num+front_num;
    if chosen_num>popu_num
        temp_num=popu_num-chosen_num+front_num;
        f(chosen_num-front_num+1:popu_num,:) = [front_i(index_of_distance(1:temp_num),:),front.*ones(temp_num,1),distance(index_of_distance(1:temp_num))];
    else
        f(chosen_num-front_num+1:chosen_num,:) = [front_i,front.*ones(front_num,1),distance];
    end
    front =  front + 1;
end
%% Crowding distance
% The crowing distance
%   Get c_sorted_based_on_front = (chromosome,order)
%   For each front Fi - c_front = c_sorted_based_on_front(1:length(Fi))
%      For each objective function - m
%         get c_sorted_based_on_objective = (c_front,m).
%         get max(m) min(m)
%         assign Inf(infinite distance) to boundary values
%         calculate the distance(m) c_front(,distance(m))
%           ? = (pre - next)/max(m) - min(m)
%      get distance sum(distance(m),2)
%
%% References
% The fast sort algorithm [1]
% [1] *Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal, and T. Meyarivan*, |A Fast
% Elitist Multiobjective Genetic Algorithm: NSGA-II|, IEEE Transactions on 
% Evolutionary Computation 6 (2002), no. 2, 182 ~ 197.
%
% the original NSGA ([2]) 
% [2] *N. Srinivas and Kalyanmoy Deb*, |Multiobjective Optimization Using 
% Nondominated Sorting in Genetic Algorithms|, Evolutionary Computation 2 
% (1994), no. 3, 221 ~ 248.
