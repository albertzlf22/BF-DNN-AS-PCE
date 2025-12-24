function y=GA_solver(x,c)
% x 个体的可优化参数值，每一行代表一个个体，每一列代表一个可优化参数
% y 计算得到的目标参数值，每一行代表一个个体，每一列代表一个目标参数
% c 用于混沌多项式代理模型的计算，对于神经网络代理模型，需重新编写该程序，可将c改为传入已训练完成的神经网络参数

%% 一维流量传感器简化模型(NN)
num=size(x,1);
y=zeros(num,2);
% tmp=bitor(x(:,16).*110+x(:,18).*10+10>120,x(:,3)*150+50+x(:,4)*8+12+x(:,1)*10+10+5>=225);
tmp=(70-40.*x(:,14)-10.*x(:,16)-15.*x(:,20)-15.*x(:,23))<0;
P=tmp.*20;    %异常几何惩罚项
y(:,1)=double(gather(extractdata(predict(c.dlnet_d,gpdl(single(x'),c.label)))))'+P;
y(:,2)=-double(gather(extractdata(predict(c.dlnet_T,gpdl(single(x'),c.label)))))'+P;
% y(:,3)=double(gather(extractdata(predict(c.dlnet_s,gpdl(single(x'),c.label)))))'+P;


function dlx = gpdl(x,labels)
dlx = dlarray(x,labels);
end

end