clear all
%dt=0.1;N=10000;Ne=8000;Ni=N-Ne;f=0.1;
D1=110;%contrast factor of connection strength
windw=0.3;% time bin
pr=1;%the part of neurons to count in

% d1=80:10:200;
% num=0;expect=0;
% for D1=d1

data0=load(['num_parameter_0_',num2str(D1),'.log']);
N=data0(1);PE=data0(2);D1=data0(3);f=data0(4);dt=data0(5);life=data0(6);
muEext=data0(7);muIext=data0(8);Tprestim_PT=data0(9);TCamp=data0(10);
Ne=round(N*PE/100);Ni=N-Ne;

% data1=load('SA_rates_pops.log');  %during SA trail
% data2=load('SA_stp_u.log');
% data3=load('SA_stp_x.log');
% data5=load('SA_currents.log');
% raster=load('SA_rasters.log');
% index_data4=find(mod(raster(:,1),10)~=0);
% data4=raster;data4(index_data4,:)=[];
% pic_all( data1,data2,data4,data3,data5,dt,Ne,f,D1);
 

% figure();
% pic_raster_et( data1,data2,data4,data3,data5,dt,Ne,f,TCamp);

data1=load(['rates_pops_0_',num2str(D1),'.log']); %after SA trail
data2=load(['stp_u_0_',num2str(D1),'.log']);
data3=load(['stp_x_0_',num2str(D1),'.log']);
data5=load(['currents_0_',num2str(D1),'.log']);   
raster=load(['rasters_0_',num2str(D1),'.log']);
index_data4=find(mod(raster(:,1),10)~=0);
data4=raster;data4(index_data4,:)=[];

figure();
pic_raster_et( data1,data2,data4,data3,data5,dt,Ne,f,D1);
% figure();
% expect=find_Power_Law(raster,windw,life,Ne,pr,f,D1 )
% pic_all( data1,data2,data4,data3,data5,dt,Ne,f,D1);

% num=num+1;
% expect(num)=find_Power_Law(raster,windw,life,Ne,pr,f,D1 );
% end
% figure();
% plot(d1./100,expect,'bo-');hold on;plot(0:0.01:2,1,'r-');hold off;grid;
% xlabel('D1','FontWeight','demi','FontSize',12);
% ylabel('\sigma','FontWeight','demi','FontSize',12);
% title('branches parameter','FontWeight','demi','FontSize',12);
% axis([0,2,0,2]);
 
