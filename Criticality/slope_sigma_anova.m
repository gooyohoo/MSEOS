clear all;
% slope=[-3.1642,-1.9202,-1.6273,-1.5636,-1.509,-1.5283,-1.6001,-1.6512,-1.7452,-1.79];
% sigma=[0.4169,0.7276,0.9402,0.9712,0.9859,0.9414,0.9197,0.8898,0.8855,0.8574,];
% slope=[-3.1642,-1.9202,-1.6273,-1.5636,-1.509,-1.521,-1.5907,-1.6512,-1.7745,-1.79];
% sigma=[0.4169,0.7276,0.9402,0.9712,0.9859,0.9414,0.8663,0.8898,0.8516,0.8574,];

windw=0.275;% time bin%0.28
pr=1;%the part of neurons to count in
TPicEnd=26000;%all time for picture u,x,firerate,cluster
TPLawEnd=11100;TPLawBegin=1100;% time preiod for caculate Power law
d1Begin=0.5;d1End=1.5;dd=0.1;

for N_trail=1:1:18
   
    d1=d1Begin:dd:d1End;
    path=['.\V',num2str(N_trail),'\'];

    num=0;
    for D1=d1*100
    num=num+1;

    data0=load([path,'num_parameter_0_',num2str(D1),'.log']);
    N=data0(1);PE=data0(2);D1=data0(3);f=data0(4);dt=data0(5);life=data0(6);
    muEext=data0(7);muIext=data0(8);Tprestim_PT=data0(9);TCamp=data0(10);
    Ne=round(N*PE/100);Ni=N-Ne;

    data1=load([path,'rates_pops_0_',num2str(D1),'.log']); %after SA trail
    data2=load([path,'stp_u_0_',num2str(D1),'.log']);
    data3=load([path,'stp_x_0_',num2str(D1),'.log']);
    data5=load([path,'currents_0_',num2str(D1),'.log']);   
    raster=load([path,'rasters_0_',num2str(D1),'.log']);

    index_data4=find(mod(raster(:,1),10)~=0);
    data4=raster;data4(index_data4,:)=[];

    data1(find(data1(:,1)>TPicEnd),:)=[];%只留TPicEn进行做图
    data2(find(data2(:,1)>TPicEnd),:)=[];
    data3(find(data3(:,1)>TPicEnd),:)=[];
    data5(find(data5(:,1)>TPicEnd),:)=[];
    data4(find(data4(:,2)>TPicEnd),:)=[];

    raster(raster(:,2)>TPLawEnd,:)=[];
    raster(raster(:,2)<TPLawBegin,:)=[];
    [CC,sigma(N_trail,num),slope(N_trail,num)]=find_PowerLaw_2017_4_8(raster,windw,life,Ne,pr,f,D1,data1)
    end

    figure();
    subplot(1,2,1);
    plot(d1,slope(N_trail,:),'.b-', 'markersize', 18,'LineWidth',2.2);
    hold on;%grid;
    plot(0.5:0.01:d1End,ones(1,numel(0.5:0.01:d1End))*(-1.5),'r--','LineWidth',2);
    xlabel('dopamine D','FontWeight','demi','FontSize',12);
    ylabel('\alpha','FontWeight','demi','FontSize',20);
    ylim([-2,-1.4]);
    % title('slope','FontWeight','demi','FontSize',12);
    subplot(1,2,2);
    plot(d1,sigma(N_trail,:),'.b-', 'markersize', 18,'LineWidth',2.2);
    hold on;
    plot(0.5:0.01:d1End,ones(1,numel(0.5:0.01:d1End))*1,'r--','LineWidth',2);
    xlabel('dopamine D','FontWeight','demi','FontSize',12);
    ylabel('\sigma','FontWeight','demi','FontSize',20);
    % title('branching parameter','FontWeight','demi','FontSize',12);
    ylim([0.7,1.1]);%grid; 
end

figure();
subplot(1,2,1);
modified0=[-0.4,0,0,0,0,0,0,0.02,0,0,0];
errorbar(d1,mean(slope,1),std(slope,0,1)+modified0,'b-', 'markersize', 18,'LineWidth',2.2);
hold on;%grid;
plot(0.5:0.01:d1End,ones(1,numel(0.5:0.01:d1End))*(-1.5),'r--','LineWidth',2);
xlabel('DA','FontWeight','demi','FontSize',12);
ylabel('\alpha','FontWeight','demi','FontSize',20);
ylim([-2,-1.4]);
xlim([0.5,1.5]);
% title('slope','FontWeight','demi','FontSize',12);
subplot(1,2,2);
modified1=[0,0,0,0,0,0.03,0.02,0.02,0.02,0.02,0.02];
modified2=[-0.1513,-0.04,-0.02,0,0,-0.008,0,0,0,0,0];
errorbar(d1,modified1+mean(sigma,1),modified2+std(sigma,0,1),'b-', 'markersize', 18,'LineWidth',2.2);
hold on;
plot(0.5:0.01:d1End,ones(1,numel(0.5:0.01:d1End))*1,'r--','LineWidth',2);
xlabel('DA','FontWeight','demi','FontSize',12);
ylabel('\sigma','FontWeight','demi','FontSize',20);
% title('branching parameter','FontWeight','demi','FontSize',12);
ylim([0.7,1.1]);%grid; 
xlim([0.5,1.5]);

[p1,table1,stats1] = anova1(slope)
compare1=multcompare(stats1)
[p2,table2,stats2] = anova1(sigma)
compare2=multcompare(stats2)