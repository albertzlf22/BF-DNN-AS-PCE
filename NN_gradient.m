function res=NN_gradient(dlnet,x,act)
% 根据网络参数求某一点的梯度（只支持一次求一个点）
% 函数自身不带归一化，需对x手动归一化!!!
if nargin==2
    act = 'relu';
end
nl=size(dlnet.Layers,1)/2;
z=cell(nl-1,1);
f=x;
for i=1:nl-1
    nn=dlnet.Layers(i*2,1).OutputSize;
    z{i}=dlnet.Layers(i*2,1).Weights*f+dlnet.Layers(i*2,1).Bias;
    if strcmp(act,'relu')
        f=max(z{i},zeros(nn,1));
    else
        f=((z{i}>0).*(1-dlnet.Layers(i*2+1,1).Scale)+dlnet.Layers(i*2+1,1).Scale).*z{i};
    end
end
% res=dlnet.Layers(nl*2,1).Weights*f+dlnet.Layers(nl*2,1).Bias;   %网络输出

% ---反向传播---
dydx=dlnet.Layers(nl*2,1).Weights;
for i=nl-1:-1:1
    if strcmp(act,'relu')
        tmp=dlnet.Layers(i*2,1).Weights.*(z{i}>0);
    else
        tmp=dlnet.Layers(i*2,1).Weights.*((z{i}>0).*(1-dlnet.Layers(i*2+1,1).Scale)+dlnet.Layers(i*2+1,1).Scale);
    end
    dydx=dydx*tmp;
end
res=dydx;
end