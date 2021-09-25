clear all
%dt=0.1;N=10000;Ne=8000;Ni=N-Ne;f=0.1;
% windw=0.4;
target=200;
d1_start=0;d1_end=2;d_d1=0.05;%contrast factor of connection strength
D1_start=d1_start*100;D1_end=d1_end*100;d_D1=d_d1*100;
D1_refer=1;


data0=load(['num_parameter_0_',num2str(D1_start),'.log']);
N=data0(1);f=data0(4);dt=data0(5);life=data0(6);PE=data0(2);
muEext=data0(7);muIext=data0(8);Tprestim_PT=data0(9);TCamp=data0(10);Tcue_PT=data0(11);
Ne=round(N*PE/100);Ni=N-Ne;

size_num=round((D1_end-D1_start)/d_D1)+1;
VC=zeros(1,size_num);
recurrent_All=zeros(1,size_num);
recurrent_E1=zeros(1,size_num);

n_nostim=round(Tprestim_PT/TCamp);
n_rest=round((life-Tprestim_PT-Tcue_PT)/TCamp);
not_count_T=round((Tprestim_PT+Tcue_PT)/TCamp);
num=0;D1=D1_start;
while D1>=D1_start&&D1<=D1_end
num=num+1;
    
raster=load(['rasters_0_',num2str(D1),'.log']);
raster(find(raster(:,2)<(Tprestim_PT+Tcue_PT)),:)=[];
T_train=raster(find(raster(:,1)==target),2);%ISI
Tisi=T_train(2:end)-T_train(1:end-1);
Pisi=tabulate(Tisi);
Pisi(find(Pisi(:,2)==0),:)=[];
VC(num)=std(Tisi)/mean(Tisi);%vari coefficient
% figure();
% plot(Pisi(:,1),Pisi(:,3),'.-');

data5=load(['currents_0_',num2str(D1),'.log']); 
data5(find(data5(:,1)<(Tprestim_PT+Tcue_PT)),:)=[];
% data5(find(data5(:,1)>=(Tprestim_PT+Tcue_PT)),:)=[];
recurrent_E1(num)=mean(data5(:,2));
recurrent_All(num)=mean(data5(:,2)+data5(:,3)+data5(:,4)+data5(:,5)+data5(:,6)+data5(:,7));


D1=D1+d_D1;
end

xx=d1_start:d_d1:d1_end;
figure();
subplot(3,1,1);
plot(xx,VC,'o-');
title('VC','FontWeight','demi','FontSize',12);
subplot(3,1,2);
plot(xx,recurrent_E1,'o-');
title('recurrent_{E1}','FontWeight','demi','FontSize',12);
subplot(3,1,3);
plot(xx,recurrent_All,'o-');
title('recurrent_{All}','FontWeight','demi','FontSize',12);

