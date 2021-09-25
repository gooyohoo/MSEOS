clear all
d1=[80:10:120,122:1:124,130:10:160];
num=0;expect=0;
for D1=d1
num=num+1;   
data0=load(['num_parameter_0_',num2str(D1),'.log']);
N=data0(1);PE=data0(2);D1=data0(3);f=data0(4);dt=data0(5);life=data0(6);
muEext=data0(7);muIext=data0(8);Tprestim_PT=data0(9);TCamp=data0(10);
Ne=round(N*PE/100);Ni=N-Ne;

data1=load(['rates_pops_0_',num2str(D1),'.log']); %after SA trail
% data2=load(['stp_u_0_',num2str(D1),'.log']);
% data3=load(['stp_x_0_',num2str(D1),'.log']);
data5=load(['currents_0_',num2str(D1),'.log']);   
raster=load(['rasters_0_',num2str(D1),'.log']);
% index_data4=find(mod(raster(:,1),10)~=0);
% data4=raster;data4(index_data4,:)=[];


fir_all=mean(data1(:,2:7));curr_all=mean(data5(:,2:7));
fir_max(num)=max(fir_all(1:5));curr_max(num)=max(curr_all(1:5));
fir_I(num)=fir_all(6);curr_I(num)=curr_all(6);
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

save('fire_train.mat','d1','fir_max','fir_I')

