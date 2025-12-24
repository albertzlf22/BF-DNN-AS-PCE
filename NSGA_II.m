function res=NSGA_II(popu_num,M,V,min_range,max_range,gen_num,solver)
% popu_num 种群数量
% M 优化目标参数数量
% V 可优化参数数量
% min_range 可优化参数范围下限，如两个可优化参数，可取最小值分别为0和1，则min_range=[0,1]
% max_range 可优化参数范围上限
% gen_num 代的总数，即总计遗传多少代后停止程序
% solver 传入一个函数，用于计算每个个体的目标参数值，参考GA_solver.m
% res 优化结果，每一行表示一个个体，前V列为优化参数值，接下来M列为优化目标值，最后两列分别为Pareto支配级别与拥挤度
init_chromosome = initialize_variables(popu_num, V, min_range, max_range);
objective = solver(init_chromosome);
popu = non_domination_sort([init_chromosome,objective], popu_num, M, V);
fig=figure;
fig_cnt=1;
scatter(objective(:,1),objective(:,2),'.');
hold on;
drawnow;
frame = getframe(fig);
im{fig_cnt} = frame2im(frame);
fig_cnt=fig_cnt+1;
% history=zeros(popu_num,M+V+2,gen_num);
for i=1:gen_num
    child_chromosome = genetic_operator(popu(:,1:V), popu(:,V+M+1:V+M+2), popu_num, 1, 20, min_range, max_range);
    chromosome = [popu(:,1:V);child_chromosome];    %精英策略
%     chromosome = child_chromosome;    %无精英策略
    objective = solver(chromosome);
    popu = non_domination_sort([chromosome,objective], popu_num, M, V);
%     history(:,:,i)=popu;
    
    if mod(i,5)==0
        scatter(popu(:,V+1),popu(:,V+2),'.');   %绘图仅能展现前两个目标参数
        drawnow;
        frame = getframe(fig);
        im{fig_cnt} = frame2im(frame);
        fig_cnt=fig_cnt+1;
    end
end
% save('history.mat','history');
res=popu;

%绘制动图
% filename = 'N_HeatTrans_Opt.gif'; 
% for i=1:fig_cnt-1
%     [A,map] = rgb2ind(im{i},256);
%     if i == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
%     end
% end

end