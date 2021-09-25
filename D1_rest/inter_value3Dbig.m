clear all
order0=0;order1=0;num=0;numy=0;x_transition=0;
% yy=[106,114,119,122,128,134,141,146,152,156,159,163];%85开始
yy=[114,119,122,128,134,141,146,152,156,159,163];%90开始

je2i=90:5:160;
figure();
for JE2I=je2i
% for JE2I=85:5:115
  storageName = strcat('.\','huhu_win_critical_Jc_Pic2_D1_EE',num2str(JE2I),'_rest','\mingw5\fire_train.mat');
  load(storageName);
  num=numel(d1);
  order1=order0+num;
%   x(order0+1:order1)=JE2I;y(order0+1:order1)=d1;z(order0+1:order1)=fir_I;  
  x(order0+1:order1)=JE2I;y(order0+1:order1)=d1;z(order0+1:order1)=fir_max;  
  plot3(x(order0+1:order1),y(order0+1:order1),z(order0+1:order1),'.-');hold on;
%   figure();plot(y(order0+1:order1),z(order0+1:order1),'.-');xlabel(num2str(JE2I));

  if JE2I<=140
  numy=numy+1;
  x_transition(numy)=order0+find(d1==yy(numy));
  end
  order0=order1;
end
%%%沿着上下相变点插值%%%
l = length(x(x_transition));
sp1 = spline(1:l,[x(x_transition);y(x_transition);z(x_transition)],1:.01:l);
sp2 = spline(1:l,[x(x_transition+1);y(x_transition+1);z(x_transition+1)],1:0.01:l);

x=[x,sp1(1,:),sp2(1,:)];y=[y,sp1(2,:),sp2(2,:)];z=[z,sp1(3,:),sp2(3,:)];


xlabel('A_{E2I}');ylabel('A_{E2E}');zlabel('fire rate[Hz]');
% title('3D fire rate transition of excitatoty');

figure();
[X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x))',linspace(min(y),max(y)),'linear');%插值
pcolor(X,Y,Z);shading interp;colorbar; %伪彩色图
xlabel('A_{E2I}');ylabel('A_{E2E}');
% title('fire rate transition of excitatoty');
axis([90,160,90,160]);

% %%%%%%%%%%%%%%此处1和后面2只能选一处%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%dopamine 22 %%%%%%%%%%%%%%
x=0:0.1:2;% 多巴胺的浓度范围
% x=[0:0.06:0.6,0.6:0.02:1.4,1.4:0.06:2];% 多巴胺的浓度范围
Cmax=1.6;Cmin=0.9;%强度因子最大和最小值
x0=1;%性能最佳的浓度
kc=0.145;%确定多巴胺调节的剧烈程度（即：对应倒U形的变化快慢）
xv=0.105;%同步调节x0e和x0i(即调xv确保（x0e+x0i）/2=x0);
x0e=x0-xv;x0i=x0+xv;%使得在x0两曲线的值等于Jei和Jee在最佳浓度的值；
Cv=Cmax-Cmin;%强度最大偏移量
Co=1-Cmin;%偏1因子

%%%%%%%%%%%%%%%%%%%%%%%%%dopamine 3 %%%%%%%%%%%%%%
% x=0:0.1:2;% 多巴胺的浓度范围
% Cmax=1.6;Cmin=1.2;%强度因子最大和最小值
% x0=1;%性能最佳的浓度
% kc=0.12;%确定多巴胺调节的剧烈程度（即：倒U形的变化快慢）
% xv=0.185;%同步调节x0e和x0i(即调xv确保（x0e+x0i）/2=x0);
% x0e=x0-xv;x0i=x0+xv;%使得在x0两曲线的值等于Jei和Jee在最佳浓度的值；
% Cv=Cmax-Cmin;%强度最大偏移量
% Co=1-Cmin;%偏1因子



xs=0;%归一化起点值到1
ce=1/(1+Cv./(1+exp((x0e-xs)./kc)));
ci=1/(1+Cv./(1+exp((x0i-xs)./kc)));

Ce=(1+Cv./(1+exp((x0e-x)./kc))).*ce-Co;%dopamine的调节
Ci=(1+Cv./(1+exp((x0i-x)./kc))).*ci-Co;
Ce_100=Ce*100;Ci_100=Ci*100;
hold on;
plot(Ci_100,Ce_100,'.-r','LineWidth',2);

%%%%%%%%%%%%%此处2和前面1只能选一处%%%%%%%%%%%%%%%%%
% figure();
% % [X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x))',linspace(min(y),max(y)),'linear');%插值
% [X,Y,Z]=griddata(x,y,z,linspace(90,160)',linspace(90,160),'linear');%插值
% surf(X,Y,Z);colorbar;%三维曲面
% xlabel('A_{E2I}','FontWeight','demi','FontSize',12);
% ylabel('A_{E2E}','FontWeight','demi','FontSize',12);
% zlabel('fire rate[Hz]','FontWeight','demi','FontSize',12);
% % title('3D fire rate transition of excitatoty');
% % axis([90,160,90,160,0,20]);

