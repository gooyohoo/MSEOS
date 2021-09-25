clear all
%dt=0.1;N=10000;Ne=8000;Ni=N-Ne;f=0.1;
% windw=0.4;
d1_start=0;d1_end=2;d_d1=0.1;%contrast factor of connection strength
D1_start=d1_start*100;D1_end=d1_end*100;d_D1=d_d1*100;
D1_refer=1;

size_num=round((D1_end-D1_start)/d_D1)+1;
u_train=zeros(1,size_num);
x_train=zeros(1,size_num);
ux_train=zeros(1,size_num);
n_T=zeros(1,size_num);

data0=load(['num_parameter_0_',num2str(D1_start),'.log']);
N=data0(1);f=data0(4);dt=data0(5);life=data0(6);PE=data0(2);
muEext=data0(7);muIext=data0(8);Tprestim_PT=data0(9);TCamp=data0(10);Tcue_PT=data0(11);
Ne=round(N*PE/100);Ni=N-Ne;

n_nostim=round(Tprestim_PT/TCamp);
n_rest=round((life-Tprestim_PT-Tcue_PT)/TCamp);
not_count_T=round((Tprestim_PT+Tcue_PT)/TCamp);
num=0;D1=D1_start;
while D1>=D1_start&&D1<=D1_end
num=num+1;
    
data1=load(['rates_pops_0_',num2str(D1),'.log']); %count n_T
data1(find(data1(:,1)<not_count_T*TCamp),:)=[];
count_t=(data1(:,2)>40&data1(:,7)>15);%确定PS的取阈值为40和15
count_t=count_t';
count=diff(count_t);
num_1=size(find(count==1),2);num_f1=size(find(count==-1),2);
if num_1-num_f1<=1
    n_T(num)=min(num_1,num_f1);
else
    n_T(num)=9999;% indicate error
end

data2=load(['stp_u_0_',num2str(D1),'.log']);
u_train(num)=mean(data2(not_count_T+1:end,2));
u_r(num)=mean(data2(1:n_nostim,2));
u_st_end(num)=data2(not_count_T,2);

data3=load(['stp_x_0_',num2str(D1),'.log']);
x_train(num)=mean(data3(not_count_T+1:end,2));
x_r(num)=mean(data3(1:n_nostim,2));
x_st_end(num)=data3(not_count_T,2);

data4=load(['rates_pops_0_',num2str(D1),'.log']);
Fe1_r(num)=mean(data4(1:n_nostim,2));
Fe1_st_end(num)=mean(data4(n_nostim:not_count_T,2));

%%%%%%%%%%%%%%%%%%%%%% figure of U（t)*X(t)%%%%%%%%%%%%%%%%%%%%%%%
% figure();
% ut=data2(1:end,2);xt=data3(1:end,2);
% ux_t=(ut.*xt);
% plot(data2(:,1),ux_t);
% line([Tprestim_PT,Tprestim_PT],[0,0.8],'color','r','line','--');
% line([Tprestim_PT+Tcue_PT,Tprestim_PT+Tcue_PT],[0,0.8],'color','r','line','--');
% title('u(t)*x(t)','FontWeight','demi','FontSize',12);
% xlabel(['t[ms]','  (D1=',num2str(D1), ')'],'FontWeight','demi','FontSize',12); 
% ylabel('Value','FontWeight','demi','FontSize',12);
% ylim([0,0.6]);grid;
  
D1=D1+d_D1;
end
ux_r=u_r.*x_r;
ux_st_end=x_st_end.*u_st_end;
ux_train=x_train+u_train;

xx=d1_start:d_d1:d1_end;

%%%%%%%%%%%%%%%%figure of U_r,X_r,UX_r,Fe1_r%%%%%%%%%
figure();
subplot(2,4,1);
plot(xx,u_r,'o-');
title('u of rest','FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
ylabel('U_r','FontWeight','demi','FontSize',12);
hold on;line([0.51,0.51],[0,0.8],'color','r','line','--');hold off;
ylim([0,0.8]);grid;
subplot(2,4,2);
plot(xx,x_r,'o-');
title('x of rest','FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
ylabel('X_r','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0.6,1.1],'color','r','line','--');hold off;
ylim([0.6,1.1]);grid;
subplot(2,4,3);
plot(xx,ux_r,'o-');
title('ux of rest','FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
ylabel('UX_r','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0,0.6],'color','r','line','--');hold off;
ylim([0,0.6]);grid;
subplot(2,4,4);
plot(xx,Fe1_r,'o-');
title('fire rate of E1 during rest','FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
ylabel('Fe1_r[Hz]','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0,7],'color','r','line','--');hold off;
ylim([0,7]);grid;
%%%%figure of U_st_end,x_st_end,UX_st_end,Fe1_st_end%%%%%%%%%%
subplot(2,4,5);
plot(xx,u_st_end,'o-');
title('u of stim end','FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
ylabel('U_stiend','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0,1],'color','r','line','--');hold off;
ylim([0,1]);grid;
subplot(2,4,6);
plot(xx,x_st_end,'o-');
title('x of stim end','FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
ylabel('X_stiend','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0,0.8],'color','r','line','--');hold off;
ylim([0,0.8]);grid;
subplot(2,4,7);
plot(xx,ux_st_end,'o-');
title('ux of stim end','FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
ylabel('UX_stiend','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0.1,0.4],'color','r','line','--');hold off;
ylim([0.1,0.4]);grid;
subplot(2,4,8);
plot(xx,Fe1_st_end,'o-');
title('fire rate of E1 when stim end','FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
ylabel('Fe1 st end[Hz]','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0,40],'color','r','line','--');hold off;
ylim([0,40]);grid;


%%%%%%%%%%%%%%%%%%%%%%%figure of mean of u,x,ux,n_T%%%%%%%%%%%%%%%%
figure();
subplot(2,2,1); 
plot(xx,x_train,'o-');
title(['-\Delta','x'],'FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
ylabel(['-\Delta','x'],'FontWeight','demi','FontSize',12);
xlim([0,2]);%ylim([0,1]); 
subplot(2,2,2);
plot(xx,u_train,'o-');
title(['\Delta','u'],'FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
ylabel(['\Delta','u'],'FontWeight','demi','FontSize',12);
xlim([0,2]);%ylim([0,1]);
subplot(2,2,3);
plot(xx,ux_train,'o-');
ylabel(['\Delta ','u+','(-\Delta','x)'],'FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
title(['\Delta ','u+','(-\Delta','x)'],'FontWeight','demi','FontSize',12);
xlim([0,2]);%ylim([0,1]);
subplot(2,2,4);
plot(xx,n_T,'o-');
title('PS number ','FontWeight','demi','FontSize',12);
xlabel(['D1','(muEext=',num2str(muEext), ')'],'FontWeight','demi','FontSize',12); 
ylabel('N_T','FontWeight','demi','FontSize',12);
xlim([0,2]);ylim([0,15]);grid;
