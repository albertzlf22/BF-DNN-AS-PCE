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
%     dt=x(i,:);
%     model.param.set('res_gap', [num2str(dt(1).*50+100),'[um]']);
%     model.param.set('hole_gap_1', [num2str(dt(2).*40+10),'[um]']);
%     model.param.set('hole_gap_2', [num2str(dt(3).*40+10),'[um]']);
%     model.param.set('x1', [num2str(dt(4).*150+50),'[um]']);
%     model.param.set('x2', [num2str(dt(5).*130+50),'[um]']);
%     model.param.set('xw1', [num2str(dt(6).*10+10),'[um]']);
%     model.param.set('xw2', [num2str(dt(7).*10+10),'[um]']);
%     model.param.set('xw3', [num2str(dt(8).*15+5),'[um]']);
%     model.param.set('xw4', [num2str(dt(9).*15+5),'[um]']);
%     model.param.set('xw5', [num2str(dt(10).*15+5),'[um]']);
%     model.param.set('xw6', [num2str(dt(11).*15+5),'[um]']);
%     model.param.set('xw7', [num2str(dt(12).*5+15),'[um]']);
%     model.param.set('xw8', [num2str(dt(13).*5+15),'[um]']);
%     model.param.set('y1', [num2str(dt(14).*40+100),'[um]']);
%     model.param.set('y2', [num2str(dt(15).*20+40),'[um]']);
%     model.param.set('yw1', [num2str(dt(16).*10+10),'[um]']);
%     model.param.set('yw2', [num2str(dt(17).*10+10),'[um]']);
%     model.param.set('yw3', [num2str(dt(18).*15+5),'[um]']);
%     model.param.set('yw4', [num2str(dt(19).*15+5),'[um]']);
%     model.param.set('yw5', [num2str(dt(20).*15+5),'[um]']);
%     model.param.set('yw6', [num2str(dt(21).*15+5),'[um]']);
%     model.param.set('yw7', [num2str(dt(22).*15+5),'[um]']);
%     model.param.set('yw8', [num2str(dt(23).*15+5),'[um]']);

%     dt=x(i,:).*0.3;
%     %opt2
%     model.param.set('res_gap', [num2str(dt(1)+150),'[um]']);
%     model.param.set('hole_gap_1', [num2str(dt(2)+20),'[um]']);
%     model.param.set('hole_gap_2', [num2str(dt(3)+20),'[um]']);
%     model.param.set('x1', [num2str(dt(4)+200),'[um]']);
%     model.param.set('x2', [num2str(dt(5)+164.94),'[um]']);
%     model.param.set('xw1', [num2str(dt(6)+16.21),'[um]']);
%     model.param.set('xw2', [num2str(dt(7)+18.38),'[um]']);
%     model.param.set('xw3', [num2str(dt(8)+20),'[um]']);
%     model.param.set('xw4', [num2str(dt(9)+20),'[um]']);
%     model.param.set('xw5', [num2str(dt(10)+19.96),'[um]']);
%     model.param.set('xw6', [num2str(dt(11)+20),'[um]']);
%     model.param.set('xw7', [num2str(dt(12)+19.97),'[um]']);
%     model.param.set('xw8', [num2str(dt(13)+19.93),'[um]']);
%     model.param.set('y1', [num2str(dt(14)+140),'[um]']);
%     model.param.set('y2', [num2str(dt(15)+35.03),'[um]']);
%     model.param.set('yw1', [num2str(dt(16)+20),'[um]']);
%     model.param.set('yw2', [num2str(dt(17)+11.805),'[um]']);
%     model.param.set('yw3', [num2str(dt(18)+10),'[um]']);
%     model.param.set('yw4', [num2str(dt(19)+19.38),'[um]']);
%     model.param.set('yw5', [num2str(dt(20)+5.79),'[um]']);
%     model.param.set('yw6', [num2str(dt(21)+8.96),'[um]']);
%     model.param.set('yw7', [num2str(dt(22)+5.01),'[um]']);
%     model.param.set('yw8', [num2str(dt(23)+5),'[um]']);
    %opt1
% %     model.param.set('res_gap', [num2str(dt(1)+150),'[um]']);
% %     model.param.set('hole_gap_1', [num2str(dt(2)+20),'[um]']);
% %     model.param.set('hole_gap_2', [num2str(dt(3)+20),'[um]']);
% %     model.param.set('x1', [num2str(dt(4)+200),'[um]']);
% %     model.param.set('x2', [num2str(dt(5)+179.65),'[um]']);
% %     model.param.set('xw1', [num2str(dt(6)+20),'[um]']);
% %     model.param.set('xw2', [num2str(dt(7)+19.30),'[um]']);
% %     model.param.set('xw3', [num2str(dt(8)+20),'[um]']);
% %     model.param.set('xw4', [num2str(dt(9)+20),'[um]']);
% %     model.param.set('xw5', [num2str(dt(10)+20),'[um]']);
% %     model.param.set('xw6', [num2str(dt(11)+20),'[um]']);
% %     model.param.set('xw7', [num2str(dt(12)+17.23),'[um]']);
% %     model.param.set('xw8', [num2str(dt(13)+19.67),'[um]']);
% %     model.param.set('y1', [num2str(dt(14)+140),'[um]']);
% %     model.param.set('y2', [num2str(dt(15)+40.1),'[um]']);
% %     model.param.set('yw1', [num2str(dt(16)+19.17),'[um]']);
% %     model.param.set('yw2', [num2str(dt(17)+11.30),'[um]']);
% %     model.param.set('yw3', [num2str(dt(18)+5.02),'[um]']);
% %     model.param.set('yw4', [num2str(dt(19)+20),'[um]']);
% %     model.param.set('yw5', [num2str(dt(20)+5),'[um]']);
% %     model.param.set('yw6', [num2str(dt(21)+9.56),'[um]']);
% %     model.param.set('yw7', [num2str(dt(22)+5.09),'[um]']);
% %     model.param.set('yw8', [num2str(dt(23)+5.21),'[um]']);
%     model.component('comp1').geom('geom1').run;
%     model.study('std1').run;
%     model.result.numerical('av1').setResult;
%     str=mphtable(model,'tbl1');
%     T1=str.data;
%     model.result.numerical('av2').setResult;
%     str=mphtable(model,'tbl3');
%     T2=str.data;
%     y(i)=T1-T2;
    %%%T型热膨胀%%%
%     dt=x(i,:);
%     w2=5+dt(1)*0.2;
%     w3=20+dt(2)*0.2;
%     th=2+dt(3)*0.06;
%     model.param.set('w2', [num2str(w2),'[um]']);
%     model.param.set('w3', [num2str(w3),'[um]']);
%     model.param.set('th', [num2str(th),'[um]']);
% 	Rt0=2e-5*200e-6/(th*w2*1e-12);
%     model.study('std1').run;
%     model.result.numerical('av1').setResult;
%     str=mphtable(model,'tbl1');
%     u=str.data;
%     model.result.numerical('av2').setResult;
%     str=mphtable(model,'tbl2');
%     v=str.data;
%     Rab0=v(1,2)/0.01e-3;
%     r1=v(2,2)*1e3/(5-2);
%     r2=v(3,2)*1e3/5;
%     u1=u(2,2);
%     u2=u(3,2);
%     Rabt=(2-u1).*(r2-r1)./(u2-u1)+r1;
% %     y_test(i)=TStruct_HeatExp(2e-6,300e-6,100e-6,5e-6,20e-6,Rabt,Rab0,Rt0,1.25e-3);
% 	y(i)=TStruct_HeatExp(2e-6,300e-6,100e-6,5e-6,20e-6,Rabt,Rab0,Rt0,1.25e-3);
end

%% Pi形残余应力
model0=HeatExp_SingleL_Epi_HF_0V;
% model0.component('comp1').mesh('mesh1').autoMeshSize(3);
err=[];
r=zeros(num,10);

for n=1:num
%     dt=xp(n,:);
    dt=x(n,:);
    w1=10+0.1*dt(1);
    w2=10+0.1*dt(2);
    w3=4+0.1*dt(3);
    L3=4+0.1*dt(4);
    w4=4+0.1*dt(5);
    L4=4+0.1*dt(6);
    L5=w4+2;
    k=28+0.5*dt(7);
    p0=7.725e-5+0.5e-5*dt(8);
    a1=1.3e-4+0.05e-4*dt(9);
    E0=138e9+1e9*dt(10);
%     w1=11.5+0.1*dt(1);
%     w2=11.5+0.1*dt(2);
%     w3=6+0.1*dt(3);
%     L3=6+0.1*dt(4);
%     w4=6+0.1*dt(5);
%     L4=6+0.1*dt(6);
%     L5=w4+1;
%     k=25+0.5*dt(7);
%     p0=7.5e-5+0.5e-5*dt(8);
%     a1=1.5e-4+0.05e-4*dt(9);
%     E0=150e9+1e9*dt(10);
    model0.param.set('w1', [num2str(w1),'[um]']);
    model0.param.set('w2', [num2str(w2),'[um]']);
    model0.param.set('w3', [num2str(w3),'[um]']);
    model0.param.set('L3', [num2str(L3),'[um]']);
    model0.param.set('w4', [num2str(w4),'[um]']);
    model0.param.set('L4', [num2str(L4),'[um]']);
    model0.param.set('L5', [num2str(L5),'[um]']);
    model0.param.set('k', num2str(k));
    model0.param.set('alpha', num2str(2e-6));
    model0.param.set('p0', num2str(p0));
    model0.param.set('a1', num2str(a1));
    model0.param.set('E0', num2str(E0));
    model0.param.set('Rs',[num2str(-10),'[MPa]']);
    try
        model0.study('std1').run;
        model0.result.numerical('av1').setResult;
        str=mphtable(model0,'tbl1');
        u0=str.data(1,2);
        model0.result.numerical('av2').setResult;
        str=mphtable(model0,'tbl2');
        R0=0.001/str.data(2,2);
        r(n,1:2)=[u0,R0];
    catch
        err=[err;n];
    end
end

clear model0;
model=HeatExp_SingleL_Epi_HF;
model.component('comp1').mesh('mesh1').autoMeshSize(4);
% model.component('comp1').physics('ec').feature('ge1').set('equation', '(aveop1(u)-ds)/0.02[um]');
for n=1:num
%     dt=xp(n,:);
    dt=x(n,:);
    u0=r(n,1);
    if u0>2
        continue;
    end
    w1=10+0.1*dt(1);
    w2=10+0.1*dt(2);
    w3=4+0.1*dt(3);
    L3=4+0.1*dt(4);
    w4=4+0.1*dt(5);
    L4=4+0.1*dt(6);
    L5=w4+2;
    k=28+0.5*dt(7);
    p0=7.725e-5+0.5e-5*dt(8);
    a1=1.3e-4+0.05e-4*dt(9);
    E0=138e9+1e9*dt(10);
%     w1=11.5+0.1*dt(1);
%     w2=11.5+0.1*dt(2);
%     w3=6+0.1*dt(3);
%     L3=6+0.1*dt(4);
%     w4=6+0.1*dt(5);
%     L4=6+0.1*dt(6);
%     L5=w4+1;
%     k=25+0.5*dt(7);
%     p0=7.5e-5+0.5e-5*dt(8);
%     a1=1.5e-4+0.05e-4*dt(9);
%     E0=150e9+1e9*dt(10);
    model.param.set('w1', [num2str(w1),'[um]']);
    model.param.set('w2', [num2str(w2),'[um]']);
    model.param.set('w3', [num2str(w3),'[um]']);
    model.param.set('L3', [num2str(L3),'[um]']);
    model.param.set('w4', [num2str(w4),'[um]']);
    model.param.set('L4', [num2str(L4),'[um]']);
    model.param.set('L5', [num2str(L5),'[um]']);
    model.param.set('k', num2str(k));
    model.param.set('alpha', num2str(2e-6));
    model.param.set('p0', num2str(p0));
    model.param.set('a1', num2str(a1));
    model.param.set('E0', num2str(E0));
    model.param.set('Rs',[num2str(-10),'[MPa]']);
    try
        model.param.set('ds',[num2str(2+u0),'[um]']);
        model.study('std1').run;
        model.result.numerical('av1').setResult;
        str=mphtable(model,'tbl1');
        u1=str.data;
        model.result.numerical('av2').setResult;
        str=mphtable(model,'tbl2');
        I1=str.data;
        model.result.numerical('gev1').setResult;
        str=mphtable(model,'tbl4');
        v1=str.data;
        model.result.numerical('av3').setResult;
        str=mphtable(model,'tbl5');
        T1=str.data;
        
        model.param.set('ds',[num2str(2-u0),'[um]']);
        model.study('std1').run;
        model.result.numerical('av1').setResult;
        str=mphtable(model,'tbl1');
        u2=str.data;
        model.result.numerical('av2').setResult;
        str=mphtable(model,'tbl2');
        I2=str.data;
        model.result.numerical('gev1').setResult;
        str=mphtable(model,'tbl4');
        v2=str.data;
        model.result.numerical('av3').setResult;
        str=mphtable(model,'tbl5');
        T2=str.data;
        
        r(n,3:10)=[u1,I1,v1,T1,u2,I2,v2,T2];
    catch
        err=[err;n];
%         flag = 0;
%         disp('Error')
%         d_in(i,:)=[];
%         d_out(i,:)=[];
    end
end


dR1=r(:,5)./r(:,4);
dR2=r(:,9)./r(:,8);
dT1=(dR1-r(:,2))./r(:,2)./1.3e-4;
dT2=(dR2-r(:,2))./r(:,2)./1.3e-4;
% y=r;
y=138e9.*2e-6.*8e-6.*(dT2-dT1)./(2.*0.41.*404e-6.*309e-6.*dT1);
% y=[r,138e9.*2e-6.*8e-6.*(dT2-dT1)./(2.*404e-6.*309e-6.*dT1)];
if ~isempty(err)
    error('There are errors')
end
end