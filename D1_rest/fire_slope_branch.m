clear all
EI_train=90:5:160;
windw=0.28;% time bin
pr=1;

d1_train1=[ 80    90   100   110   114   115   116   120   130   140   150   160];%90
d1_train2=[80    90   100   110   116   117   118   119   120   121   122   123   130   140   150   160];%95
d1_train3=[80    90   100   110   120   122   123   124   130   140   150   160];%100
d1_train4=[80    90   100   110   120   125   126   127   128   129   130   140   150   160];%105
d1_train5=[ 80    90   100   110   120   130   134   135   136   140   150   160];%110
d1_train6=[80    90   100   110   120   130   140   141   142   143   150   160];%115
d1_train7=[80    90   100   110   120   130   140   146   147   150   160];%120
d1_train8=[ 80    90   100   110   120   130   140   151   152   153   154   160];%125
d1_train9= [80    90   100   110   120   130   140   150   156   157   158   160];%130
d1_train10=[ 80    90   100   110   120   130   140   150   158   159   162   163];%135
d1_train11=[80    90   100   110   120   130   140   150   160   162   163   164];%140
d1_train12= [80    90   100   110   120   130   140   150   160];%145
d1_train13=[ 80    90   100   110   120   130   140   150   160];%150
d1_train14=[ 80    90   100   110   120   130   140   150   160];%155
d1_train15=[ 80    90   100   110   120   130   140   150   160];%160
numb=0;
for EI=EI_train
numb=numb+1;

storageName = strcat('.\huhu_win_critical_Jc_Pic2_D1_EE',num2str(EI),'_rest','\mingw5');
eval(['d1=d1_train',num2str(numb),';']);
disp(EI)
disp(d1)

num=0;fir_max=0;curr_max=0;fir_I=0;curr_I=0;expect=0;slope=0;
for D1=d1
num=num+1;   
data0=load([storageName,'\num_parameter_0_',num2str(D1),'.log']);
N=data0(1);PE=data0(2);D1=data0(3);f=data0(4);dt=data0(5);life=data0(6);
muEext=data0(7);muIext=data0(8);Tprestim_PT=data0(9);TCamp=data0(10);
Ne=round(N*PE/100);Ni=N-Ne;

data1=load([storageName,'\rates_pops_0_',num2str(D1),'.log']); %after SA trail
% data2=load(['stp_u_0_',num2str(D1),'.log']);
% data3=load(['stp_x_0_',num2str(D1),'.log']);
data5=load([storageName,'\currents_0_',num2str(D1),'.log']);   
raster=load([storageName,'\rasters_0_',num2str(D1),'.log']);
% index_data4=find(mod(raster(:,1),10)~=0);
% data4=raster;data4(index_data4,:)=[];


fir_all=mean(data1(:,2:7));curr_all=mean(data5(:,2:7));

fir_max(num)=max(fir_all(1:5));curr_max(num)=max(curr_all(1:5));
fir_I(num)=fir_all(6);curr_I(num)=curr_all(6);
%figure(011111111);
[CC,expect(num),slope(num)]=find_Power_Law(raster,windw,life,Ne,pr,f,D1,data1 );

end

figure();
subplot(2,2,1);
plot(d1./100,fir_max,'ro-');hold on;
xlabel('D1','FontWeight','demi','FontSize',12);
ylabel('Hz','FontWeight','demi','FontSize',12);
title('fire rate(max) of E','FontWeight','demi','FontSize',12);
subplot(2,2,2);
plot(d1./100,fir_I,'bo-');hold on;
xlabel('D1','FontWeight','demi','FontSize',12);
ylabel('Hz','FontWeight','demi','FontSize',12);
title('fire rate of I','FontWeight','demi','FontSize',12);
subplot(2,2,3);
plot(d1./100,curr_max,'ro-');hold on;
xlabel('D1','FontWeight','demi','FontSize',12);
ylabel('mA','FontWeight','demi','FontSize',12);
title('recurrent(max) of E','FontWeight','demi','FontSize',12);
subplot(2,2,4);
plot(d1./100,curr_I,'bo-');hold on;
xlabel('D1','FontWeight','demi','FontSize',12);
ylabel('mA','FontWeight','demi','FontSize',12);
title('recurrent of I','FontWeight','demi','FontSize',12);

figure();
plot(d1./100,expect,'bo-');hold on;plot(0:0.01:2,1,'r-');hold off;grid;
xlabel('D1','FontWeight','demi','FontSize',12);
ylabel('\sigma','FontWeight','demi','FontSize',12);
title('branches parameter','FontWeight','demi','FontSize',12);
% axis([0,2,0,2]);
figure();
plot(d1./100,slope,'bo-');hold on;plot(0:0.01:2,1,'r-');hold off;grid;
xlabel('D1','FontWeight','demi','FontSize',12);
ylabel('slope','FontWeight','demi','FontSize',12);
title('power law','FontWeight','demi','FontSize',12);

save([storageName,'\slope_train.mat'],'d1','slope');
%save([storageName,'\fire_train.mat'],'d1','fir_max','fir_I')
save([storageName,'\branches_train.mat'],'d1','expect');
end

