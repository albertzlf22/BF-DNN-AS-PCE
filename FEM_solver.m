function y=FEM_solver(model,x)

num = size(x,1);
% dim = 23;
% if nargin == 2
%     vec = eye(dim);
% end
% % x = [randn(num,dim-1),x]/vec;
% % x = [zeros(num,dim-1),x]/vec;
% % x = [x(:,1),zeros(num,dim-2),x(:,2)]/vec;
% x = [x(:,1),randn(num,dim-2),x(:,2)]/vec;
y = zeros(num,1);
for i=1:num
    %%%流量%%%
    dt=x(i,:).*0.3;

    %opt1
    model.param.set('res_gap', [num2str(dt(1)+150),'[um]']);
    model.param.set('hole_gap_1', [num2str(dt(2)+20),'[um]']);
    model.param.set('hole_gap_2', [num2str(dt(3)+20),'[um]']);
    model.param.set('x1', [num2str(dt(4)+200),'[um]']);
    model.param.set('x2', [num2str(dt(5)+179.65),'[um]']);
    model.param.set('xw1', [num2str(dt(6)+20),'[um]']);
    model.param.set('xw2', [num2str(dt(7)+19.30),'[um]']);
    model.param.set('xw3', [num2str(dt(8)+20),'[um]']);
    model.param.set('xw4', [num2str(dt(9)+20),'[um]']);
    model.param.set('xw5', [num2str(dt(10)+20),'[um]']);
    model.param.set('xw6', [num2str(dt(11)+20),'[um]']);
    model.param.set('xw7', [num2str(dt(12)+17.23),'[um]']);
    model.param.set('xw8', [num2str(dt(13)+19.67),'[um]']);
    model.param.set('y1', [num2str(dt(14)+140),'[um]']);
    model.param.set('y2', [num2str(dt(15)+40.1),'[um]']);
    model.param.set('yw1', [num2str(dt(16)+19.17),'[um]']);
    model.param.set('yw2', [num2str(dt(17)+11.30),'[um]']);
    model.param.set('yw3', [num2str(dt(18)+5.02),'[um]']);
    model.param.set('yw4', [num2str(dt(19)+20),'[um]']);
    model.param.set('yw5', [num2str(dt(20)+5),'[um]']);
    model.param.set('yw6', [num2str(dt(21)+9.56),'[um]']);
    model.param.set('yw7', [num2str(dt(22)+5.09),'[um]']);
    model.param.set('yw8', [num2str(dt(23)+5.21),'[um]']);
    model.component('comp1').geom('geom1').run;
    model.study('std1').run;
    model.result.numerical('av1').setResult;
    str=mphtable(model,'tbl1');
    T1=str.data;
    model.result.numerical('av2').setResult;
    str=mphtable(model,'tbl3');
    T2=str.data;
    y(i)=T1-T2;
end

end