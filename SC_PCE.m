function c=SC_PCE(d,order,solver,typ,algor)
% 
if nargin==4
    if d<=3
        algor='TP';
    else
        algor='SG';
    end
end

if strcmp(typ,'Hermite')
    typ_str = 'GQN';
else
    typ_str = 'GQU';
end

if strcmp(algor,'TP')
    [x_tmp,w_tmp] = nwspgr(typ_str, 1, order+1, 1);
    ind=(1:order+1)';
    for i=2:d
        ind=[kron((1:order+1)',ones((order+1)^(i-1),1)),kron(ones(order+1,1),ind)];
    end
    x = x_tmp(ind);
    w = prod(w_tmp(ind),2);
else
    [x,w] = nwspgr(typ_str, d, order+1, 1);
end

if strcmp(typ,'Legendre')
    x = x.*2-1;
%     w = w.*2;
end

[F,~,~] = basis_table(typ,order);

y = solver(x);
t_d = total_degree(d,order);
total_basis = size(t_d,1);
p_num = length(w);
A = ones(total_basis,p_num);
for i=1:total_basis
    for j=1:d
        A(i,:)=A(i,:).*F{t_d(i,j)+1}(x(:,j))';
    end
end
c = A*(w.*y);

if strcmp(typ,'Legendre')
    p = 1./(2.*(0:order)'+1);
else
    p = factorial(0:order)';
end
tmp = prod(p(t_d+1),2);
c = c./tmp;

end