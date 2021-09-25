clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%dopamine 2 %%%%%%%%%%%%%%
x=0:0.1:2;% 多巴胺的浓度范围
% x=[0:0.06:0.6,0.6:0.02:1.4,1.4:0.06:2];% 多巴胺的浓度范围
Cmax=1.6;Cmin=0.9;%强度因子最大和最小值
x0=1;%性能最佳的浓度
kc=0.15;%确定多巴胺调节的剧烈程度（即：对应倒U形的变化快慢）
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
% % kc=0.15;%确定多巴胺调节的剧烈程度（即：对应倒U形的变化快慢）
% x0e=x0-xv;x0i=x0+xv;%使得在x0两曲线的值等于Jei和Jee在最佳浓度的值；
% Cv=Cmax-Cmin;%强度最大偏移量
% Co=1-Cmin;%偏1因子


xs=0;%归一化起点值到1
ce=1/(1+Cv./(1+exp((x0e-xs)./kc)));
ci=1/(1+Cv./(1+exp((x0i-xs)./kc)));

xstim=0:0.01:2;
Ce=(1+Cv./(1+exp((x0e-x)./kc))).*ce-Co;%dopamine的调节
Ci=(1+Cv./(1+exp((x0i-x)./kc))).*ci-Co;
Ce_100=Ce*100;Ci_100=Ci*100;
Ce2=(1+Cv./(1+exp((x0e-xstim)./kc))).*ce-Co;%dopamine的调节
Ci2=(1+Cv./(1+exp((x0i-xstim)./kc))).*ci-Co;
Ce2_100=Ce2*100;Ci2_100=Ci2*100;

figure();
% % subplot(1,2,1)
% plot(x,Ce_100,'r-',x,Ci_100,'b--','LineWidth',2.5);%grid;
% xlabel('Dopamine','FontSize',12);ylabel('A','FontSize',12);
% plot(xstim,Ce2_100,xstim,Ci2_100,xstim,(Ce2_100-Ci2_100),'LineWidth',2.5);
% % xlabel('Ci');ylabel('Ce');
% legend('A_{E2E}','A_{E2I}','FontSize',12);
% 
DeltaA=Ce2_100-Ci2_100;
[hAx,hLine1,hLine2] = plotyy(xstim,[Ci2_100;Ce2_100],xstim,DeltaA);
xlabel('DA');
ylabel(hAx(1),'A [%]') ;
ylim(hAx(1),[Cmin*100,Cmax*100]);
ylim(hAx(2),[0,max(DeltaA)]);
ylabel(hAx(2),'\DeltaA [%]');
legend('A_{EI}','A_{EE}','\Delta A','FontSize',12);
set(hLine1, 'linestyle', '-','LineWidth',2.5 ); 
set(hLine2, 'linestyle', '--','LineWidth',2.5 ); 

% subplot(1,2,2)
figure();
% x=80:5:140;y=[103,109,114,117,123,128,135,139,145,148,155,162,163];%上界
% % x=80:5:140;y=[103,108,114,119,122,128,134,141,146,152,156,162,163];%下届
x=90:5:140;y=[114,119,122,128,134,141,146,152,156,162,163];%下届

l = length(x);
sp1 = spline(1:l,[x;y],1:.1:l);%插值

plot(sp1(1,:),sp1(2,:),'r','LineWidth',2.5,'MarkerSize',17.5);hold on;plot(80:0.1:160,80:0.1:160,'k--','LineWidth',1.5);%grid;
% plot(x,y,'r.',sp1(1,:),sp1(2,:),'r','LineWidth',1.8,'MarkerSize',17.5);hold on;plot(80:0.1:160,80:0.1:160,'k--','LineWidth',1.5);grid;
xlabel('A_{E2I}','FontSize',12);ylabel('A_{E2E}','FontSize',12);
axis([90,160,90,160]);
hold on;
plot(Ci_100,Ce_100,'.-b');
legend('transition','reference','dopamine','FontSize',12)
grid;