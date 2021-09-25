clear all
order0=0;order1=0;num=0;numy=0;x_transition=0;
% yy=[106,114,119,122,128,134,141,146,152,156,159,163];%85��ʼ
yy=[114,119,122,128,134,141,146,152,156,159,163];%90��ʼ

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
%%%�������������ֵ%%%
l = length(x(x_transition));
sp1 = spline(1:l,[x(x_transition);y(x_transition);z(x_transition)],1:.01:l);
sp2 = spline(1:l,[x(x_transition+1);y(x_transition+1);z(x_transition+1)],1:0.01:l);

x=[x,sp1(1,:),sp2(1,:)];y=[y,sp1(2,:),sp2(2,:)];z=[z,sp1(3,:),sp2(3,:)];


xlabel('A_{E2I}');ylabel('A_{E2E}');zlabel('fire rate[Hz]');
% title('3D fire rate transition of excitatoty');

figure();
[X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x))',linspace(min(y),max(y)),'linear');%��ֵ
pcolor(X,Y,Z);shading interp;colorbar; %α��ɫͼ
xlabel('A_{E2I}');ylabel('A_{E2E}');
% title('fire rate transition of excitatoty');
axis([90,160,90,160]);

% %%%%%%%%%%%%%%�˴�1�ͺ���2ֻ��ѡһ��%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%dopamine 22 %%%%%%%%%%%%%%
x=0:0.1:2;% ��Ͱ���Ũ�ȷ�Χ
% x=[0:0.06:0.6,0.6:0.02:1.4,1.4:0.06:2];% ��Ͱ���Ũ�ȷ�Χ
Cmax=1.6;Cmin=0.9;%ǿ������������Сֵ
x0=1;%������ѵ�Ũ��
kc=0.145;%ȷ����Ͱ����ڵľ��ҳ̶ȣ�������Ӧ��U�εı仯������
xv=0.105;%ͬ������x0e��x0i(����xvȷ����x0e+x0i��/2=x0);
x0e=x0-xv;x0i=x0+xv;%ʹ����x0�����ߵ�ֵ����Jei��Jee�����Ũ�ȵ�ֵ��
Cv=Cmax-Cmin;%ǿ�����ƫ����
Co=1-Cmin;%ƫ1����

%%%%%%%%%%%%%%%%%%%%%%%%%dopamine 3 %%%%%%%%%%%%%%
% x=0:0.1:2;% ��Ͱ���Ũ�ȷ�Χ
% Cmax=1.6;Cmin=1.2;%ǿ������������Сֵ
% x0=1;%������ѵ�Ũ��
% kc=0.12;%ȷ����Ͱ����ڵľ��ҳ̶ȣ�������U�εı仯������
% xv=0.185;%ͬ������x0e��x0i(����xvȷ����x0e+x0i��/2=x0);
% x0e=x0-xv;x0i=x0+xv;%ʹ����x0�����ߵ�ֵ����Jei��Jee�����Ũ�ȵ�ֵ��
% Cv=Cmax-Cmin;%ǿ�����ƫ����
% Co=1-Cmin;%ƫ1����



xs=0;%��һ�����ֵ��1
ce=1/(1+Cv./(1+exp((x0e-xs)./kc)));
ci=1/(1+Cv./(1+exp((x0i-xs)./kc)));

Ce=(1+Cv./(1+exp((x0e-x)./kc))).*ce-Co;%dopamine�ĵ���
Ci=(1+Cv./(1+exp((x0i-x)./kc))).*ci-Co;
Ce_100=Ce*100;Ci_100=Ci*100;
hold on;
plot(Ci_100,Ce_100,'.-r','LineWidth',2);

%%%%%%%%%%%%%�˴�2��ǰ��1ֻ��ѡһ��%%%%%%%%%%%%%%%%%
% figure();
% % [X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x))',linspace(min(y),max(y)),'linear');%��ֵ
% [X,Y,Z]=griddata(x,y,z,linspace(90,160)',linspace(90,160),'linear');%��ֵ
% surf(X,Y,Z);colorbar;%��ά����
% xlabel('A_{E2I}','FontWeight','demi','FontSize',12);
% ylabel('A_{E2E}','FontWeight','demi','FontSize',12);
% zlabel('fire rate[Hz]','FontWeight','demi','FontSize',12);
% % title('3D fire rate transition of excitatoty');
% % axis([90,160,90,160,0,20]);

