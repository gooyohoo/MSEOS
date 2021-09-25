%论文中用的就是EIbalance3
clear all
DTrain=0:25:200;
% DTrain=100;
DopamineLong=numel(DTrain);
TargetTrain1=8:8:800;%E1
TargetTrain2=1608:8:2400;%E3
TargetTrain3=4040:40:8000;%En
TargetTrain4=8020:20:10000;%I

TargetTrain=TargetTrain1;%inital

TargLong=numel(TargetTrain);
TargetNum=4;
Cepre=zeros(TargetNum,DopamineLong,TargLong);
Cedur=zeros(TargetNum,DopamineLong,TargLong);
Ceaft=zeros(TargetNum,DopamineLong,TargLong);
  
num2=0;
for D1=DTrain
% D1=100;%contrast factor of connection strength
num2=num2+1;

data0=load(['num_parameter_0_',num2str(D1),'.log']);
N=data0(1);PE=data0(2);f=data0(4);dt=data0(5);life=data0(6);
muEext=data0(7);muIext=data0(8);Tprestim_PT=data0(9);TCamp=data0(10);Tcue_PT=data0(11);
Ne=round(N*PE/100);Ni=N-Ne;

data6=load(['I_Exc_0_',num2str(D1),'.log']); %after SA trail
data7=load(['I_Inh_0_',num2str(D1),'.log']);
data8=load(['I_Sum_0_',num2str(D1),'.log']);
data9=load(['I_ExtBg_0_',num2str(D1),'.log']);
data11=load(['I_Leaky_0_',num2str(D1),'.log']);
plook=1;TCamp1=10;
Lmax=life/TCamp1;
Tpre=1:Tprestim_PT/TCamp1;Tdur=Tprestim_PT/TCamp1+1:(Tprestim_PT+Tcue_PT)/TCamp1;
Taft=(Tprestim_PT+Tcue_PT)/TCamp1+1:life/TCamp1;

% figure();
% xx=(1:life/plook)*TCamp1;xx1=(1:life)*TCamp1;
% subplot(2,1,1);
% plot(xx,data6(1:end/plook,target)+data9(1:end/plook,target),'r');
% title(['I_{Exc+Ext} of cell ',num2str(target)],'FontWeight','demi','FontSize',12);
% xlabel('T[ms]','FontWeight','demi','FontSize',12); 
% grid;
% subplot(2,1,2);
% plot(xx,-data7(1:end/plook,target),'b');
% title(['- I_{Inh} of cell ',num2str(target)],'FontWeight','demi','FontSize',12);
% xlabel('T[ms]','FontWeight','demi','FontSize',12); 
% grid;

for num1=1:1:4
    eval(['Target=TargetTrain',num2str(num1),';']);
  

    num3=0;
    for target=Target
        num3=num3+1;

        Data11=data6(:,target)%+data9(:,target);
        Data12=data7(:,target);
%         Data12=data7(:,target)+data11(:,target);
        
        Xpre=Data11(Tpre,:);Xdur=Data11(Tdur,:);Xaft=Data11(Taft,:);
        Ypre=Data12(Tpre,:);Ydur=Data12(Tdur,:);Yaft=Data12(Taft,:);

        th1=-10;th2=-10;%去峰
        Xpre(find(Ypre<th1))=[];Xdur(find(Ydur<th1))=[];Xaft(find(Yaft<th2))=[];
        Ypre(find(Ypre<th1))=[];Ydur(find(Ydur<th1))=[];Yaft(find(Yaft<th2))=[];

        Cepre(num1,num2,num3)=min(min(corrcoef(Xpre,Ypre)));
        Cedur(num1,num2,num3)=min(min(corrcoef(Xdur,Ydur)));
        Ceaft(num1,num2,num3)=min(min(corrcoef(Xaft,Yaft)));    
    end
end
disp([num2str(D1),'->',num2str(DTrain(end))]);
end

MeanCepre=zeros(TargetNum,DopamineLong);MeanCedur=zeros(TargetNum,DopamineLong);MeanCeaft=zeros(TargetNum,DopamineLong);
StdCepre=zeros(TargetNum,DopamineLong);StdCedur=zeros(TargetNum,DopamineLong);StdCeaft=zeros(TargetNum,DopamineLong);
figure();
for num=1:1:4
eval(['TargetTrain=TargetTrain',num2str(num),';']);   

% figure();
% MeanCepre(num,:)=mean(Cepre(num,:,:),3);MeanCedur(num,:)=mean(Cedur(num,:,:),3);MeanCeaft(num,:)=mean(Ceaft(num,:,:),3);
% StdCepre(num,:)=std(Cepre(num,:,:),0,3);StdCedur(num,:)=std(Cedur(num,:,:),0,3);StdCeaft(num,:)=std(Ceaft(num,:,:),0,3);
% subplot(1,3,1);
% errorbar(DTrain*0.01,MeanCepre(num,:),StdCepre(num,:),'-o');
% title('Mean[Cpre]','FontWeight','demi','FontSize',12);xlim([-0.5,2.5]);
% subplot(1,3,2);
% errorbar(DTrain*0.01,MeanCedur(num,:),StdCedur(num,:),'-o');
% xlabel(['T[ms]',' pop ',num2str(TargetTrain(1))],'FontWeight','demi','FontSize',12);
% title('Mean[Cdur]','FontWeight','demi','FontSize',12);xlim([-0.5,2.5]);
% subplot(1,3,3);
% errorbar(DTrain*0.01,MeanCeaft(num,:),StdCeaft(num,:),'-o');
% title('Mean[Caft]','FontWeight','demi','FontSize',12);xlim([-0.5,2.5]);

MeanCepre(num,:)=mean(Cepre(num,:,:),3);MeanCedur(num,:)=mean(Cedur(num,:,:),3);MeanCeaft(num,:)=mean(Ceaft(num,:,:),3);
StdCepre(num,:)=std(Cepre(num,:,:),0,3);StdCedur(num,:)=std(Cedur(num,:,:),0,3);StdCeaft(num,:)=std(Ceaft(num,:,:),0,3);
subplot(1,4,num);
errorbar(DTrain*0.01,-MeanCedur(num,:),StdCedur(num,:),'-o');
title('Mean[Cdur]','FontWeight','demi','FontSize',12);xlim([-0.5,2.5]);ylim([0,1])
xlabel([' pop ',num2str(TargetTrain(1))],'FontWeight','demi','FontSize',12);

% figure();
% subplot(1,3,1);
% plot(DTrain*0.01,squeeze(Cepre(num,:,:))','.');
% title('Cpre','FontWeight','demi','FontSize',12);
% subplot(1,3,2);
% plot(DTrain*0.01,squeeze(Cedur(num,:,:))','.');
% title('Cdur','FontWeight','demi','FontSize',12);
% subplot(1,3,3);
% plot(DTrain*0.01,squeeze(Ceaft(num,:,:))','.');
% title('Caft','FontWeight','demi','FontSize',12);
% xlabel(['T[ms]',' pop ',num2str(TargetTrain(1))],'FontWeight','demi','FontSize',12);
end


%%
%Anova test
for popN=1:1:4
   [p,c,s]=anova1((squeeze(Cedur(popN,:,:)))');%去掉；可以得到boxplot 看到方差逐渐缩小；倒U形曲线的中间方差小，也是他们趋同的结果；
   PValue(popN)=p;  
   compare(:,:,popN)=multcompare(s);
end
 

%after
% figure();
% for num=1:4
% MeanCepre(num,:)=mean(Cepre(num,:,:),3);MeanCedur(num,:)=mean(Cedur(num,:,:),3);MeanCeaft(num,:)=mean(Ceaft(num,:,:),3);
% StdCepre(num,:)=std(Cepre(num,:,:),0,3);StdCedur(num,:)=std(Cedur(num,:,:),0,3);StdCeaft(num,:)=std(Ceaft(num,:,:),0,3);
%   
% subplot(1,4,num);
% errorbar(DTrain*0.01,-MeanCeaft(num,:),StdCeaft(num,:),'-o');
% title('Mean[Cdur]','FontWeight','demi','FontSize',12);xlim([-0.5,2.5]);ylim([0,1])
% xlabel([' pop ',num2str(TargetTrain(1))],'FontWeight','demi','FontSize',12);
% end

%pre
% figure();
% for num=1:4
% MeanCepre(num,:)=mean(Cepre(num,:,:),3);MeanCedur(num,:)=mean(Cedur(num,:,:),3);MeanCeaft(num,:)=mean(Ceaft(num,:,:),3);
% StdCepre(num,:)=std(Cepre(num,:,:),0,3);StdCedur(num,:)=std(Cedur(num,:,:),0,3);StdCeaft(num,:)=std(Ceaft(num,:,:),0,3);
%   
% subplot(1,4,num);
% errorbar(DTrain*0.01,-MeanCepre(num,:),StdCepre(num,:),'-o');
% title('Mean[Cdur]','FontWeight','demi','FontSize',12);xlim([-0.5,2.5]);ylim([0,1])
% xlabel([' pop ',num2str(TargetTrain(1))],'FontWeight','demi','FontSize',12);
% end



