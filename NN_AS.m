function [v,e]=NN_AS(dlnet,x,dim_ind,weight,act)
if nargin==4
    act='relu';
end
n=size(x,2);
d=length(dim_ind);
C=zeros(d);
for i=1:n
    grad=NN_gradient(dlnet,x(:,i),act);
    tmp=grad(dim_ind)./weight;
    C=C+tmp'*tmp;
end
C=C./n;
[v,e]=eig(C);

% grad=NN_gradient(dlnet,x);
% tmp=grad(dim_ind)./weight;
% C=tmp'*tmp;
% [v,e]=eig(C);
end