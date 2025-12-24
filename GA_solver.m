function y=GA_solver(x,c)
% x 个体的可优化参数值，每一行代表一个个体，每一列代表一个可优化参数
% y 计算得到的目标参数值，每一行代表一个个体，每一列代表一个目标参数
% c 用于混沌多项式代理模型的计算，对于神经网络代理模型，需重新编写该程序，可将c改为传入已训练完成的神经网络参数
%% Test Function
% y=zeros(size(x,1),2);
% y(:,1)=2.*((x(:,1)-0.5).^2+(x(:,2)-0.5).^2);
% y(:,2)=2.*((x(:,1)-0.25).^2+(x(:,2)-0.25).^2);
%% CTE test structure PCE model
% [P,~,~]=basis_table('Legendre');
% t_d=total_degree(10,4);
% total_basis=size(t_d,1);
% num=size(x,1);
% y=zeros(num,2);
% r1=zeros(num,5);
% r2=zeros(num,5);
% for n=1:5
%     d1=[x(:,1:5),(-1.5+n*0.5).*ones(num,1),zeros(num,2),x(:,6:7)];
%     Mt1=ones(num,total_basis);
%     d2=[x(:,1:5),zeros(num,2),(-1.5+n*0.5).*ones(num,1),x(:,6:7)];
%     Mt2=ones(num,total_basis);
%     for i=1:num
%         for j=1:total_basis
%             for m=1:10
%                 Mt1(i,j)=Mt1(i,j)*P{t_d(j,m)+1}(d1(i,m));
%                 Mt2(i,j)=Mt2(i,j)*P{t_d(j,m)+1}(d2(i,m));
%             end
%         end
%     end
%     %直接对处理值回归
% %     r1(:,n)=(Mt1*c-25.*ones(num,1)).^2;
% %     r2(:,n)=(Mt2*c-(10+5.*n).*ones(num,1)).^2;
%     %对测量值回归
%     Vh_r=Mt1*c(:,1);
%     Vh_r0=Mt1*c(:,2);
%     Vt_r=Mt1*c(:,3);
%     Vh_t=Mt1*c(:,4);
%     Vh_t0=Mt1*c(:,5);
%     Vt_t=Mt1*c(:,6);
% %     temp1(:,n)=[Vh_r,Vh_r0,Vt_r,Vh_t,Vh_t0,Vt_t]';
%     r1(:,n)=(Heat_Trans_Pi(Vh_r,Vh_r0,Vt_r,Vh_t,Vh_t0,Vt_t,(d1+1)./2)-25.*ones(num,1)).^2;
% %     temp1(n)=Heat_Trans_Pi(Vh_r,Vh_r0,Vt_r,Vh_t,Vh_t0,Vt_t,(d1+1)./2);
%     Vh_r=Mt2*c(:,1);
%     Vh_r0=Mt2*c(:,2);
%     Vt_r=Mt2*c(:,3);
%     Vh_t=Mt2*c(:,4);
%     Vh_t0=Mt2*c(:,5);
%     Vt_t=Mt2*c(:,6);
% %     temp2(:,n)=[Vh_r,Vh_r0,Vt_r,Vh_t,Vh_t0,Vt_t]';
%     r2(:,n)=(Heat_Trans_Pi(Vh_r,Vh_r0,Vt_r,Vh_t,Vh_t0,Vt_t,(d2+1)./2)-(10+5.*n).*ones(num,1)).^2;
% %     temp2(n)=Heat_Trans_Pi(Vh_r,Vh_r0,Vt_r,Vh_t,Vh_t0,Vt_t,(d2+1)./2);
% end
% y(:,1)=mean(r1,2);
% y(:,2)=mean(r2,2);
% % y=[r1,r2];
% % y=[temp1,temp2];
% % y=temp2;
%% New CTE test structure PCE model
% [P,~,~]=basis_table('Legendre');
% t_d=total_degree(8,4);
% total_basis=size(t_d,1);
% num=size(x,1);
% y=zeros(num,2);
% r1=zeros(num,5);
% r2=zeros(num,5);
% for n=1:5
%     d1=[x(:,1:6),(-1.5+n*0.5).*ones(num,1),-0.875.*ones(num,1)];
%     Mt1=ones(num,total_basis);
%     d2=[x(:,1:6),zeros(num,1),(-1.5+n*0.5).*ones(num,1)];
%     Mt2=ones(num,total_basis);
%     for i=1:num
%         for j=1:total_basis
%             for m=1:8
%                 Mt1(i,j)=Mt1(i,j)*P{t_d(j,m)+1}(d1(i,m));
%                 Mt2(i,j)=Mt2(i,j)*P{t_d(j,m)+1}(d2(i,m));
%             end
%         end
%     end
%     k1=Mt1*c(:,1);
%     b1=Mt1*c(:,2);
%     d_in=(d1+1)./2;
%     for i=1:num
%         A=2e-6*(5e-6+d_in(i,1)*10e-6);
%         l_zh=20e-6+d_in(i,2).*30e-6;
%         l_zv=150e-6+d_in(i,3).*100e-6;
%         l=6*l_zh+4*l_zv;
%         [r1(i,n),~,~]=NStruc_HeatTrans(k1(i),b1(i),0.7e-3+0.1e-3.*n,1.1e-5,-5e-4,A,l);
%     end
%     k2=Mt2*c(:,1);
%     b2=Mt2*c(:,2);
%     d_in=(d2+1)./2;
%     for i=1:num
%         A=2e-6*(5e-6+d_in(i,1)*10e-6);
%         l_zh=20e-6+d_in(i,2).*30e-6;
%         l_zv=150e-6+d_in(i,3).*100e-6;
%         l=6*l_zh+4*l_zv;
%         [r2(i,n),~,~]=NStruc_HeatTrans(k2(i),b2(i),1e-3,1.1e-5,-5e-4,A,l);
%     end
% %     r1(:,n)=(Mt1*c-28.3.*ones(num,1)).^2;
% %     r2(:,n)=(Mt2*c-28.3.*ones(num,1)).^2;
% %     r1(:,n)=(r1(:,n)-28.3.*ones(num,1)).^2;
% %     r2(:,n)=(r2(:,n)-28.3.*ones(num,1)).^2;
% %     temp1(n)=Mt1*c;
% end
% % y(:,1)=mean(r1,2);
% % y(:,2)=mean(r2,2);
% % y=temp1;
% y=r1;
%% New CTE Test Structure Theory Analysis
% num=size(x,1);
% d=[(x+1)./2,ones(num,3),0.5.*ones(num,1),1.*ones(num,1)];
% [k,b,~,~]=NStruc_HeatTrans_Theory(d);
% y=zeros(num,2);
% for i=1:num
%     A=2e-6*(5e-6+d(i,1)*10e-6);
%     l_zh=20e-6+d(i,2).*30e-6;
%     l_zv=150e-6+d(i,3).*100e-6;
%     l=6*l_zh+4*l_zv;
%     temp1=NStruc_HeatTrans(0.99*k(i),b(i),1e-3,1.1e-5,-5e-4,A,l);
%     temp2=NStruc_HeatTrans(1.01*k(i),b(i),1e-3,1.1e-5,-5e-4,A,l);
%     y(i,1)=abs(temp1-temp2);
%     temp1=NStruc_HeatTrans(k(i),0.99*b(i),1e-3,1.1e-5,-5e-4,A,l);
%     temp2=NStruc_HeatTrans(k(i),1.01*b(i),1e-3,1.1e-5,-5e-4,A,l);
%     y(i,2)=abs(temp1-temp2);
% end
%% Wire Analysis
% [P,~,~]=basis_table('Legendre');
% t_d=total_degree(5,4);
% total_basis=size(t_d,1);
% num=size(x,1);
% y=zeros(num,4);
% r1=zeros(num,5);
% r2=zeros(num,5);
% r3=zeros(num,5);
% r4=zeros(num,5);
% for n=1:5
%     d1=[x,[(-1.5+n*0.5),-7/8].*ones(num,2)];
%     Mt1=ones(num,total_basis);
%     d2=[x,[0,(-1.5+n*0.5)].*ones(num,2)];
%     Mt2=ones(num,total_basis);
%     for i=1:num
%         for j=1:total_basis
%             for m=1:5
%                 Mt1(i,j)=Mt1(i,j)*P{t_d(j,m)+1}(d1(i,m));
%                 Mt2(i,j)=Mt2(i,j)*P{t_d(j,m)+1}(d2(i,m));
%             end
%         end
%     end
%     [ktk,btk,~,~]=NStruc_HeatTrans_Theory((d1+1)./2);
%     [kth,bth,~,~]=NStruc_HeatTrans_Theory((d2+1)./2);
%     r1(:,n)=(Mt1*c(:,1)-ktk).^2;
%     r2(:,n)=(Mt1*c(:,2)-btk).^2;
%     r3(:,n)=(Mt2*c(:,1)-kth).^2;
%     r4(:,n)=(Mt2*c(:,2)-bth).^2;
% end
% y(:,1)=mean(r1,2);
% y(:,2)=mean(r2,2);
% y(:,3)=mean(r3,2);
% y(:,4)=mean(r4,2);
%% New CTE test structure PCE model 0303 HeatTrans_N7.xlsx
% [P,~,~]=basis_table('Legendre');
% t_d=total_degree(9,4);
% total_basis=size(t_d,1);
% num=size(x,1);
% y=zeros(num,2);
% r1=zeros(num,5);
% r2=zeros(num,5);
% for n=1:5
%     d1=[x,(-1.5+n*0.5).*ones(num,1),-(2/3).*ones(num,1)];
%     Mt1=ones(num,total_basis);
%     d2=[x,zeros(num,1),(-1.5+n*0.5).*ones(num,1)];
%     Mt2=ones(num,total_basis);
%     for i=1:num
%         for j=1:total_basis
%             for m=1:9
%                 Mt1(i,j)=Mt1(i,j)*P{t_d(j,m)+1}(d1(i,m));
%                 Mt2(i,j)=Mt2(i,j)*P{t_d(j,m)+1}(d2(i,m));
%             end
%         end
%     end
%     k1=Mt1*c(:,1);
%     b1=Mt1*c(:,2);
%     d_in=(d1+1)./2;
%     for i=1:num
%         A=2e-6*(5e-6+d_in(i,1)*5e-6);
%         l_zh=20e-6+d_in(i,2).*30e-6;
%         l_zv=150e-6+d_in(i,3).*100e-6;
%         l=6*l_zh+4*l_zv;
%         [r1(i,n),~,~]=NStruc_HeatTrans(k1(i),b1(i),1e-3,1.1e-5,-5e-4,A,l);
%     end
%     k2=Mt2*c(:,1);
%     b2=Mt2*c(:,2);
%     d_in=(d2+1)./2;
%     for i=1:num
%         A=2e-6*(5e-6+d_in(i,1)*5e-6);
%         l_zh=20e-6+d_in(i,2).*30e-6;
%         l_zv=150e-6+d_in(i,3).*100e-6;
%         l=6*l_zh+4*l_zv;
%         [r2(i,n),~,~]=NStruc_HeatTrans(k2(i),b2(i),1e-3,1.1e-5,-5e-4,A,l);
%     end
% %     r1(:,n)=(Mt1*c-28.3.*ones(num,1)).^2;
% %     r2(:,n)=(Mt2*c-28.3.*ones(num,1)).^2;
%     r1(:,n)=(r1(:,n)-(10+n*5).*ones(num,1)).^2;
%     r2(:,n)=(r2(:,n)-25.*ones(num,1)).^2;
% %     temp1(n)=Mt1*c;
% end
% y(:,1)=mean(r1,2);
% y(:,2)=mean(r2,2);
% % y=temp1;
% % y=r1;
%% New CTE test structure PCE model 0420 HeatTrans_D1.xlsx
% [P,~,~]=basis_table('Legendre');
% t_d=total_degree(12,4);
% total_basis=size(t_d,1);
% num=size(x,1);
% y=zeros(num,2);
% r1=zeros(num,5);
% r2=zeros(num,5);
% for n=1:5
%     d1=[x(:,1:7),(-1.5+n*0.5).*ones(num,1),-(2/3).*ones(num,1),x(:,8:10)];
%     Mt1=ones(num,total_basis);
%     d2=[x(:,1:7),zeros(num,1),(-1.5+n*0.5).*ones(num,1),x(:,8:10)];
%     Mt2=ones(num,total_basis);
%     for i=1:num
%         for j=1:total_basis
%             for m=1:12
%                 Mt1(i,j)=Mt1(i,j)*P{t_d(j,m)+1}(d1(i,m));
%                 Mt2(i,j)=Mt2(i,j)*P{t_d(j,m)+1}(d2(i,m));
%             end
%         end
%     end
%     k1=Mt1*c(:,1);
%     b1=Mt1*c(:,2);
%     d_in=(d1+1)./2;
%     for i=1:num
%         A=2e-6*(5e-6+d_in(i,1)*5e-6);
%         l_zh=20e-6+d_in(i,2).*30e-6;
%         l_zv=150e-6+d_in(i,3).*100e-6;
%         l=6*l_zh+4*l_zv;
%         [r1(i,n),~,~]=NStruc_HeatTrans(k1(i),b1(i),1e-3,1.1e-5,-5e-4,A,l);
%     end
%     k2=Mt2*c(:,1);
%     b2=Mt2*c(:,2);
%     d_in=(d2+1)./2;
%     for i=1:num
%         A=2e-6*(5e-6+d_in(i,1)*5e-6);
%         l_zh=20e-6+d_in(i,2).*30e-6;
%         l_zv=150e-6+d_in(i,3).*100e-6;
%         l=6*l_zh+4*l_zv;
%         [r2(i,n),~,~]=NStruc_HeatTrans(k2(i),b2(i),1e-3,1.1e-5,-5e-4,A,l);
%     end
% %     r1(:,n)=(Mt1*c-28.3.*ones(num,1)).^2;
% %     r2(:,n)=(Mt2*c-28.3.*ones(num,1)).^2;
%     r1(:,n)=(r1(:,n)-(10+n*5).*ones(num,1)).^2;
%     r2(:,n)=(r2(:,n)-25.*ones(num,1)).^2;
% %     temp1(n)=Mt1*c;
% end
% y(:,1)=mean(r1,2);
% y(:,2)=mean(r2,2);
% % y=temp1;
% % y=r1;
%% Circle Ris
% num=size(x,1);
% y=zeros(num,2);
% r1=zeros(num,36);
% r2=zeros(num,36);
% act_v=10.^(-5+(0:0.2:1).*2)./5e-6;
% for p=1:6
%     dr=normrnd(0,0.6,2000,2);
%     dr(:,1)=dr(:,1)./100;
%     dr(:,2)=dr(:,2)./125;
% 	for q=1:3
%         para=[0.8*q-1.8,0.4*p-1.4];
%         temp=UQ_CircleRis(x,dr,c(:,1),para,act_v(q));
%         r1(:,p*6+q-6)=temp(:,1);
%         r2(:,p*6+q-6)=temp(:,2);
% 	end
% 	for q=1:3
%         para=[0.8*q-1.4,0.4*p-1.4];
%         temp=UQ_CircleRis(x,dr,c(:,2),para,act_v(q+3));
%         r1(:,p*6+q-3)=temp(:,1);
%         r2(:,p*6+q-3)=temp(:,2);
% 	end
% end
% y(:,1)=mean(r1,2);
% y(:,2)=mean(r2,2);
% % dr=normrnd(0,0.6,2000,2);
% % dr(:,1)=dr(:,1)./100;
% % dr(:,2)=dr(:,2)./125;
% % y=UQ_CircleRis(x,dr,c(:,2),[0.2,0],10^(-5+0.8*2)./5e-6);
%% 改进的热导率测试结构(NN)
% c为包含三个网络信息（k,b,T(用于计算灵敏度））的struct
% num=size(x,1);
% y=zeros(num,2);
% % r1=zeros(num,5);
% % r2=zeros(num,5);
% x=(x+1)./2;
% % dt1=[kron(x(:,1:7),ones(25,1)),repmat([kron((0:0.25:1)',ones(5,1)),repmat((0.25:0.125:0.75)',5,1)],num,1)];
% % kt1=double(gather(extractdata(predict(c.dlnet_k,gpdl(single(dt1'),c.label)))));
% % bt1=double(gather(extractdata(predict(c.dlnet_b,gpdl(single(dt1'),c.label)))));
% % dt2=[kron(x(:,1:7),ones(5,1)),repmat([(0:0.25:1)',0.5.*ones(5,1)],num,1)];
% % Tt2=double(gather(extractdata(predict(c.dlnet_T,gpdl(single(dt2'),c.label)))));
% % Tt2=Tt2.*c.sT+c.mT;
% dt1=[kron(x(:,1:7),ones(5,1)),repmat([(0:0.25:1)',0.5.*ones(5,1)],num,1)];
% dt2=[kron(x(:,1:7),ones(5,1)),repmat([0.5.*ones(5,1),(0.25:0.125:0.75)'],num,1)];
% kt1=double(gather(extractdata(predict(c.dlnet_k,gpdl(single(dt1'),c.label)))));
% kt2=double(gather(extractdata(predict(c.dlnet_k,gpdl(single(dt2'),c.label)))));
% bt1=double(gather(extractdata(predict(c.dlnet_b,gpdl(single(dt1'),c.label)))));
% bt2=double(gather(extractdata(predict(c.dlnet_b,gpdl(single(dt2'),c.label)))));
% for i=1:num
%     tmp1=0;
%     tmp2=0;
%     A=5e-6*(5+x(i,1)*5)*1e-6;
%     L=(6*(20+x(i,2)*30)+4*(150+x(i,3)*100))*1e-6;
%     for j=1:5
% %     for j=1:25
%         k1=kt1(i*5+j-5).*c.sk+c.mk;
%         b1=bt1(i*5+j-5).*c.sb+c.mb;
%         [res,~,~]=NStruc_HeatTrans(k1,b1,1e-3,1e-4,-5e-3,A,L);
% %         k_h=120+floor(j/5)*5;
% %         tmp1=tmp1+abs((res-k_h)/k_h)*100;
%         tmp1=tmp1+abs((res-(115+5*j))./(115+5*j))*100;
%         k2=kt2(i*5+j-5).*c.sk+c.mk;
%         b2=bt2(i*5+j-5).*c.sb+c.mb;
%         [res,~,~]=NStruc_HeatTrans(k2,b2,1e-3,1e-4,-5e-3,A,L);
%         tmp2=tmp2+abs((res-130)./130)*100;
%     end
%     y(i,1)=tmp1./5;
%     y(i,2)=tmp2./5;
% %     y(i,1)=tmp1./25;
% %     tmp2=Tt2(5*i-4:5*i);
% %     y(i,2)=-mean(abs(tmp2(2:5)-tmp2(1:4))./5);
% end
%% 改进的热导率测试结构(gpc)
% num=size(x,1);
% y=zeros(num,2);
% [P,~,~]=basis_table('Legendre');
% t_d=total_degree(9,4);
% total_basis=size(t_d,1);
% dt1=[kron(x(:,1:7),ones(5,1)),repmat([(0:0.25:1)',0.5.*ones(5,1)],num,1).*2-1];
% dt2=[kron(x(:,1:7),ones(5,1)),repmat([0.5.*ones(5,1),(0.25:0.125:0.75)'],num,1).*2-1];
% % dt1=[kron(x(:,1:7),ones(25,1)),repmat([kron((0:0.25:1)',ones(5,1)),repmat((0.25:0.125:0.75)',5,1)],num,1)];
% % dt2=[kron(x(:,1:7),ones(5,1)),repmat([(0:0.25:1)',0.5.*ones(5,1)],num,1)];
% % Mt1=ones(num*25,total_basis);
% Mt1=ones(num*5,total_basis);
% Mt2=ones(num*5,total_basis);
% for j=1:total_basis
%     for m=1:9
%         Mt1(:,j)=Mt1(:,j).*P{t_d(j,m)+1}(dt1(:,m));
%         Mt2(:,j)=Mt2(:,j).*P{t_d(j,m)+1}(dt2(:,m));
%     end
% end
% kt1=Mt1*c(:,1);
% bt1=Mt1*c(:,2);
% kt2=Mt2*c(:,1);
% bt2=Mt2*c(:,2);
% % Tt2=Mt2*c(:,3);
% for i=1:num
%     tmp1=0;
%     tmp2=0;
%     A=5e-6*(5+((x(i,1)+1)/2)*5)*1e-6;
%     L=(6*(20+((x(i,2)+1)/2)*30)+4*(150+((x(i,3)+1)/2)*100))*1e-6;
%     for j=1:5
% %     for j=1:25
%         k1=kt1(i*5+j-5);
%         b1=bt1(i*5+j-5);
%         [res,~,~]=NStruc_HeatTrans(k1,b1,1e-3,1e-4,-5e-3,A,L);
% %         k_h=120+floor((j-1)/5)*5;
% %         tmp1=tmp1+abs((res-k_h)/k_h)*100;
%         tmp1=tmp1+abs((res-(115+5*j))./(115+5*j))*100;
%         k2=kt2(i*5+j-5);
%         b2=bt2(i*5+j-5);
%         [res,~,~]=NStruc_HeatTrans(k2,b2,1e-3,1e-4,-5e-3,A,L);
%         tmp2=tmp2+abs((res-130)./130)*100;
%     end
%     y(i,1)=tmp1./5;
%     y(i,2)=tmp2./5;
% %     y(i,1)=tmp1./25;
% %     tmp2=Tt2(5*i-4:5*i);
% %     y(i,2)=-mean(abs(tmp2(2:5)-tmp2(1:4))./5).*200;
% end
%% T形热膨胀(NN)
% num=size(x,1);
% y=zeros(num,3);
% dt1=[repmat(normrnd(0,0.15,2000,3)./[5,10,15],num,1)+kron(x,ones(2000,1)),repmat(rand(2000,3),num,1),ones(num*2000,1).*0.6];
% yt1=double(gather(extractdata(predict(c.dlnet,gpdl(single(dt1'),c.label)))));
% dt2=[kron(x,ones(5,1)),ones(num*5,3).*[0.5,0.5,0.25],repmat((0:0.25:1)',num,1)];
% yt2=double(gather(extractdata(predict(c.dlnet,gpdl(single(dt2'),c.label)))));
% yt3=double(gather(extractdata(predict(c.dlnet_rt,gpdl(single(dt2'),c.label)))));
% for i=1:num
%     tmp=yt1(i*2000-1999:i*2000).*c.s+c.m;
%     y(i,1)=std(tmp);
%     tmp=yt2(i*5-4:i*5).*c.s+c.m;
%     y(i,2)=mean(abs(tmp-(2e-6:0.125e-6:2.5e-6)));
%     tmp=yt3(i*5-4:i*5).*c.s_rt+c.m_rt;
%     y(i,3)=-mean(tmp(1:4)-tmp(2:5));
% end
%% 一维流量传感器简化模型(NN)
% num=size(x,1);
% y=zeros(num,2);
% % tmp=bitor(x(:,16).*110+x(:,18).*10+10>120,x(:,3)*150+50+x(:,4)*8+12+x(:,1)*10+10+5>=225);
% tmp=(70-40.*x(:,14)-10.*x(:,16)-15.*x(:,20)-15.*x(:,23))<0;
% P=tmp.*20;    %异常几何惩罚项
% y(:,1)=double(gather(extractdata(predict(c.dlnet_d,gpdl(single(x'),c.label)))))'+P;
% y(:,2)=-double(gather(extractdata(predict(c.dlnet_T,gpdl(single(x'),c.label)))))'+P;
% % y(:,3)=double(gather(extractdata(predict(c.dlnet_s,gpdl(single(x'),c.label)))))'+P;
%% 二维流量传感器简化模型(NN)
% num=size(x,1);
% y=zeros(num,3);
% % tmp=bitor(x(:,16).*110+x(:,18).*10+10>120,x(:,3)*150+50+x(:,4)*8+12+x(:,1)*10+10+5>=225);
% tmp=x(:,16).*110+x(:,18).*10+10>120;
% P=tmp.*20;    %异常几何惩罚项
% y=double(gather(extractdata(predict(c.dlnet,gpdl(single(x'),c.label)))))'+P;
% % % y(:,1)=double(gather(extractdata(predict(c.dlnet_s,gpdl(single(x'),c.label)))))';
% % y(:,1)=double(gather(extractdata(predict(c.dlnet_w,gpdl(single(x'),c.label)))))'+P;
% % y(:,2)=-double(gather(extractdata(predict(c.dlnet_T,gpdl(single(x'),c.label)))))'+P;
% % y(:,3)=double(gather(extractdata(predict(c.dlnet_s,gpdl(single(x'),c.label)))))'+P;
%% 二维流量传感器简化模型(NN)
% num=size(x,1);
% y=zeros(num,3);
% % y=zeros(num,2);
% y(:,1)=double(gather(extractdata(predict(c.dlnet_s,gpdl(single(x'),c.label)))))';
% y(:,2)=-double(gather(extractdata(predict(c.dlnet_T1,gpdl(single(x'),c.label)))))';
% tmp=double(gather(extractdata(predict(c.dlnet_T2,gpdl(single(x'),c.label)))))';
% dT1=-y(:,2).*c.s_T1+c.m_T1;
% dT2=tmp.*c.s_T2+c.m_T2;
% y(:,3)=dT2-dT1;
%% 无刻蚀孔流量传感器简化模型(NN)
% num=size(x,1);
% y=zeros(num,2);
% x(:,8)=floor(x(:,8).*6)./5;     %离散化
% x(:,9)=floor(x(:,9).*11)./10;   %离散化
% 
% w_h=x(:,3).*12+8;
% w_hc=x(:,5).*(w_h-8)+8;
% n_h_max=floor(390./(2.*w_h+2.*w_hc));
% n_h_min=ceil(n_h_max./2);
% n_h_max=(n_h_max-2)./10;
% n_h_min=(n_h_min-2)./10;
% w_m=x(:,2).*12+8;
% w_mc=x(:,4).*(w_m-8)+8;
% gap=x(:,6).*80+20;
% L_hr=x(:,7).*30+50;
% n_m_max=floor((325-gap-L_hr)./(2.*w_m+2.*w_mc));
% n_m_min=max(ceil(n_m_max./2),2);
% n_m_max=(n_m_max-2)./5;
% n_m_min=(n_m_min-2)./5;
% tmp=bitor(x(:,8)<n_m_min,x(:,8)>n_m_max);
% tmp=bitor(tmp,bitor(x(:,9)<n_h_min,x(:,9)>n_h_max));
% P=tmp.*20;    %异常几何惩罚项
% % y(:,1)=double(gather(extractdata(predict(c.dlnet_p,gpdl(single(x'),c.label)))))'+P;
% y(:,1)=double(gather(extractdata(predict(c.dlnet_T,gpdl(single(x'),c.label)))))'+P;
% y(:,2)=-double(gather(extractdata(predict(c.dlnet_r,gpdl(single(x'),c.label)))))'+P;
%% 无刻蚀孔流量传感器加热电阻简化模型(NN)
% num=size(x,1);
% % y=zeros(num,2);
% x(:,10)=min(floor(x(:,10).*6)./5,1);     %离散化
% P=zeros(num,1);
% for i=1:num
%     dt=x(i,:);
% 
%     a_w=dt(4)*6-3;
%     tmp_min=min([(1+a_w)^2,(1-a_w)^2,0]);
%     tmp_max=max([(1+a_w)^2,(1-a_w)^2,0]);
%     c_w=(dt(5)*2-1)*(15/(tmp_max-tmp_min));
%     if c_w>0
%         b_w=dt(6).*(15-tmp_max*c_w+tmp_min*c_w)+(5-tmp_min*c_w);
%     else
%         b_w=dt(6).*(15-tmp_min*c_w+tmp_max*c_w)+(5-tmp_max*c_w);
%     end
% 
%     a_g=dt(7)*6-3;
%     tmp_min=min([(1+a_g)^2,(1-a_g)^2,0]);
%     tmp_max=max([(1+a_g)^2,(1-a_g)^2,0]);
%     c_g=(dt(8)*2-1)*(15/(tmp_max-tmp_min));
%     if c_g>0
%         b_g=dt(9).*(15-tmp_max*c_g+tmp_min*c_g)+(5-tmp_min*c_g);
%     else
%         b_g=dt(9).*(15-tmp_min*c_g+tmp_max*c_g)+(5-tmp_max*c_g);
%     end
% 
%     n_w=dt(10).*5+5;
%     W=c_w.*((-1:2/(n_w-1):1)-a_w).^2+b_w;
%     G=c_g.*((-1:2/(n_w-1):1)-a_g).^2+b_g;
%     if sum(W)+sum(G)-G(1)/2>200
%         P(i)=20;
%     end
% end
% % y=-double(gather(extractdata(predict(c.dlnet,gpdl(single(x'),c.label)))))'+[P,P];
% 
% y(:,1)=-double(gather(extractdata(predict(c.dlnet_mi,gpdl(single(x'),c.label)))))'+P;
% y(:,2)=-double(gather(extractdata(predict(c.dlnet_mn,gpdl(single(x'),c.label)))))'+P;
%% Glass 1D v1
% num=size(x,1);
% % y=zeros(num,3);
% % P=((x(:,3)*130+70)/2+x(:,1)*130+70+x(:,4)*110+70>(x(:,6)*200+600)/2-20).*20;
% % y(:,1)=double(gather(extractdata(predict(c.dlnet_P,gpdl(single(x'),c.label)))))'+P;
% % y(:,2:3)=double(gather(extractdata(predict(c.dlnet_T,gpdl(single(x'),c.label)))))'+[P,-P];
% % y(:,3)=-y(:,3);
% y=zeros(num,2);
% % P=((x(:,3)*130+70)/2+x(:,1)*130+70+x(:,4)*110+70>(x(:,6)*200+600)/2-20).*20;    %v1
% P=((x(:,3)*50+70)/2+(x(:,1)*50+70)*2+x(:,4)*50+50+x(:,5)*50+50>(x(:,7)*200+700)/2-20).*20;    %v2
% y(:,1)=double(gather(extractdata(predict(c.dlnet_P,gpdl(single(x'),c.label)))))'+P;
% y(:,2)=-double(gather(extractdata(predict(c.dlnet_dT,gpdl(single(x'),c.label)))))'+P;
%% 无刻蚀孔流量传感器加热电阻简化模型(NN) v1_1
% num=size(x,1);
% y=zeros(num,2);
% x(:,14)=min(floor(x(:,14).*8)./7,1);     %离散化
% P=zeros(num,1);
% for i=1:num
%     dt=x(i,:);
% 
%     [a_L,b_L,c_L,~]=RandQuadFuncGen(20,100,3,-3,dt(5:7));
%     [a_w,b_w,c_w,~]=RandQuadFuncGen(5,20,3,-3,dt(8:10));
%     [a_g,b_g,c_g,~]=RandQuadFuncGen(5,20,3,-3,dt(11:13));
% 
%     n_w=dt(14).*7+5;
%     n_L=ceil(n_w/2);
%     W=c_w.*((-1:2/(n_w-1):1)-a_w).^2+b_w;
%     G=c_g.*((-1:2/(n_w-1):1)-a_g).^2+b_g;
%     L=c_L.*((-1:2/(n_L-1):1)-a_L).^2+b_L;
%     heat_h=2*(sum(W)+sum(G))-G(1);
%     heat_w=max(L)*2+10;
% 
%     if heat_h>500 || heat_w/2+dt(2)*110+70+dt(1)*130+70>(dt(4)*200+600)/2-20
%         P(i)=20;
%     end
% end
% % y=-double(gather(extractdata(predict(c.dlnet,gpdl(single(x'),c.label)))))'+[P,P];
% 
% y(:,1)=-double(gather(extractdata(predict(c.dlnet_T,gpdl(single(x'),c.label)))))'+P;
% y(:,2)=double(gather(extractdata(predict(c.dlnet_P,gpdl(single(x'),c.label)))))'+P;
%% 无刻蚀孔流量传感器加热电阻简化模型(NN)-固定膜尺寸 v1_2
% num=size(x,1);
% y=zeros(num,2);
% x(:,12)=min(floor(x(:,12).*8)./7,1);     %离散化
% P=zeros(num,1);
% for i=1:num
%     dt=x(i,:);
% 
%     [a_L,b_L,c_L,~]=RandQuadFuncGen(20,80,-3,3,dt(3:5));
%     [a_w,b_w,c_w,~]=RandQuadFuncGen(5,20,-3,3,dt(6:8));
%     [a_g,b_g,c_g,~]=RandQuadFuncGen(5,20,-3,3,dt(9:11));
% 
%     n_w=dt(12).*7+5;
%     n_L=ceil(n_w/2);
%     W=c_w.*((-1:2/(n_w-1):1)-a_w).^2+b_w;
%     G=c_g.*((-1:2/(n_w-1):1)-a_g).^2+b_g;
%     L=c_L.*((-1:2/(n_L-1):1)-a_L).^2+b_L;
%     heat_h=2*(sum(W)+sum(G))-G(1);
%     heat_w=max(L)*2+10;
% 
%     if heat_h>500 || heat_w/2+dt(2)*80+50+dt(1)*80+50>(600)/2-20
%         P(i)=20;
%     end
% end
% % y=-double(gather(extractdata(predict(c.dlnet,gpdl(single(x'),c.label)))))'+[P,P];
% 
% y(:,1)=-double(gather(extractdata(predict(c.dlnet_T,gpdl(single(x'),c.label)))))'+P;
% y(:,2)=double(gather(extractdata(predict(c.dlnet_P,gpdl(single(x'),c.label)))))'+P;
%% 无刻蚀孔流量传感器加热电阻简化模型(NN) v2_1
% num=size(x,1);
% y=zeros(num,2);
% x(:,15)=min(floor(x(:,15).*8)./7,1);     %离散化
% P=zeros(num,1);
% for i=1:num
%     dt=x(i,:);
% 
%     [a_L,b_L,c_L,~]=RandQuadFuncGen(20,55,3,-3,dt(6:8));
%     [a_w,b_w,c_w,~]=RandQuadFuncGen(5,20,3,-3,dt(9:11));
%     [a_g,b_g,c_g,~]=RandQuadFuncGen(5,20,3,-3,dt(12:14));
% 
%     n_w=dt(15).*7+5;
%     n_L=ceil(n_w/2);
%     W=c_w.*((-1:2/(n_w-1):1)-a_w).^2+b_w;
%     G=c_g.*((-1:2/(n_w-1):1)-a_g).^2+b_g;
%     L=c_L.*((-1:2/(n_L-1):1)-a_L).^2+b_L;
%     heat_h=2*(sum(W)+sum(G))-G(1);
%     heat_w=max(L)*2+10;
% 
%     if heat_h>500 || heat_w/2+(dt(1)*50+70)*2+dt(2)*50+50+dt(5)*50+50>(dt(4)*200+700)/2-20
%         P(i)=20;
%     end
% end
% % y=-double(gather(extractdata(predict(c.dlnet,gpdl(single(x'),c.label)))))'+[P,P];
% 
% y(:,1)=-double(gather(extractdata(predict(c.dlnet_T,gpdl(single(x'),c.label)))))'+P;
% y(:,2)=double(gather(extractdata(predict(c.dlnet_P,gpdl(single(x'),c.label)))))'+P;
%% 谐振固支梁 GAN_AS
% dim=17;
% 
% num=size(x,1);
% y=zeros(num,3);
% dlZ=dlarray([x,zeros(num,100-dim)]/c.v,'BSC');
% dlX=GenModel(dlZ,c.paraGenerator,c.GenState,c.hyperparameters,false);
% y(:,1)=-double(extractdata(c.f_E(dlX)))';
% y(:,2)=-double(extractdata(c.f_S(dlX)))';
% y(:,3)=-double(extractdata(c.f_P(dlX)))';
%% 谐振固支梁 GAN_AS
% dim=8;
% 
% num=size(x,1);
% y=zeros(num,2);
% dlZ=dlarray([x,zeros(num,100-dim)]/c.v,'BSC');
% dlX=GenModel(dlZ,c.paraGenerator,c.GenState,c.hyperparameters,false);
% dlX(:,1,:)=dlX(:,1,:)+(-1:0.08:1)';
% y(:,1)=-double(extractdata(c.f_o(dlX)))';
% y(:,2)=double(extractdata(c.f_T(dlX)))';
%% V形梁参数化优化
num=size(x,1);
y=zeros(num,2);
y(:,1)=-double(gather(extractdata(predict(c.dlnet_P,gpdl(single([x,[0.125,0.5].*ones(num,2)]'),c.label)))))';
y(:,2)=double(gather(extractdata(predict(c.dlnet_T,gpdl(single([x,[0.125,0.5].*ones(num,2)]'),c.label)))))';
end

function res=UQ_CircleRis(x,dr,ct,para,act_v)
    [P,~,~]=basis_table('Legendre');
    t_d=total_degree(4,4);
    total_basis=size(t_d,1);
    sam_num=size(x,1);
    rnd_num=size(dr,1);
    res=zeros(sam_num,2);
    for n=1:sam_num
        d=[x(n,:)+dr,para.*ones(rnd_num,2)];
        Mt=ones(rnd_num,total_basis);
%         for i=1:rnd_num
        for j=1:total_basis
            for m=1:4
                Mt(:,j)=Mt(:,j).*P{t_d(j,m)+1}(d(:,m));
            end
        end
%         end
        r_mc=abs(Mt*ct-act_v)./act_v;
        res(n,1)=mean(r_mc);
        res(n,2)=std(r_mc);
    end
end

function dlx = gpdl(x,labels)
dlx = dlarray(x,labels);
end