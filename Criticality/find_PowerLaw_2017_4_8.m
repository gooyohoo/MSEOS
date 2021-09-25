function [CC,expect,slope ] = find_PowerLaw_2017_4_8( raster,windw,life,Ne,pr,f,D1,data1 )
%FIND_POWERLAW_2017_4_8 Summary of this function goes here
%   Detailed explanation goes here
p=5;
SMA=1000;SMI=20;SNUM1=30;SNUM2=30;NLDw=3;NLUp=180;

fir_all=mean(data1(:,2:6));
maxP=find(fir_all==max(fir_all));

% raster(find(raster(:,1)>Ne),:)=[];%只取兴奋性E细胞
% raster(find(raster(:,1)>=Ne*p*f),:)=[];%只取兴奋性E细胞的选择性细胞
% raster(find(raster(:,1)>Ne*f),:)=[];%只取兴奋性E1细胞的选择性细胞
raster(find(raster(:,1)>=Ne*maxP*f),:)=[];%只取兴奋性Emax细胞
raster(find(raster(:,1)<Ne*(maxP-1)*f),:)=[];%只取兴奋性Emax细胞

index_r=find(mod(raster(:,1),round(1/pr))~=0);%从目标集中取pr*100*%部分观察
raster(index_r,:)=[];

CC=0;
% for t=1:1:floor(life/windw) %%测试是否一个domain里面有某个神经元出现多次
% 
% AA=raster;
% AA(find(AA(:,2)>windw*(t)),:)=[];
% AA(find(AA(:,2)<windw*(t-1)),:)=[];
% BB=tabulate(AA(:,1));
% BB(find(BB(:,3)==0),:)=[];
% CC(t)=numel(find(BB(:,2)>1));
% end

raster_t=raster(:,2);
% raster_t=raster(1:round(end/10),2);
life=floor(life/windw)*windw;
central_window=(windw/2):windw:life;
fire_num=hist(raster_t,central_window);%求每个窗口中频数
fire_num=[0,fire_num,0];%确保首尾为0

index_0=find(fire_num==0);   
index_inter=index_0(:,2:size(index_0,2))-index_0(:,1:size(index_0,2)-1);%求0元素间间隔
index_inter(size(index_inter,2))=index_inter(size(index_inter,2))+1;%最后一个簇把最后0含上
P_cell=mat2cell(fire_num,1,index_inter);

P_nc=cell(size(P_cell));
for i=1:size(index_inter,2) %求和
   P_nc(i)={sum(cat(2,P_cell{1,i},0))};   
end
   P_num=cell2mat(P_nc);
   P_num(P_num==0)=[];%去0
 
%   P_num( find(P_num>500))=[];
% dlmwrite(['avalanches_',num2str(D1),'.txt'],P_num);  
  
   
 m=0;n_max=max(P_num);  
for k=1:size(index_0,2)-1
    if fire_num(index_0(k)+1)~=0
         m=m+1;
%           expect(m)=(fire_num(index_0(k)+2)/fire_num(index_0(k)+1))*(n_max-1)/(n_max-fire_num(index_0(k)+1));
%           expect(m)=round(fire_num(index_0(k)+2)/fire_num(index_0(k)+1))*(n_max-1)/(n_max-fire_num(index_0(k)+1));
         expect(m)=(fire_num(index_0(k)+2)/fire_num(index_0(k)+1));
         
    end
end
expect=mean(expect);

% subplot(1,3,1); %no average
% NUM=tabulate(P_num');
% NUM(find(NUM(:,3)==0),:)=[];
% loglog(NUM(:,1),NUM(:,3)./100,'gd-', 'markersize', 5); 
% axis([0,SMA,10^-5,1.1]);
% hold on;
bb=logspace(0,log10(SMA),SNUM1);
% loglog(bb(1:end-1), bb(1:end-1).^-1.5,'r--');grid;
% hold on;



% subplot(1,3,2); %average in log axis ALL
% bb=logspace(0,log10(SMA),SNUM1);
% del=floor(bb(2:end)-bb(1:end-1))+1;
% [n,x]=histc(P_num',bb);
% cc=sqrt(bb(2:end).*bb(1:end-1)); loglog(cc, n(1:end-1)'./del/size(P_num',1),'b.-', 'markersize', 5); 
% hold on; loglog(bb(1:end-1), bb(1:end-1).^-1.5,'r-');
% sumP2=sum( n(1:end-1)'./del/size(P_num',1))
% axis([0,SMA,10^-5,1.1]);
% hold on;

% subplot(1,3,3); %average in log axis PART
NUM=tabulate(P_num');
NUM(find(NUM(:,3)==0),:)=[];
NUM(find(NUM(:,1)>SMI),:)=[];%front
bb=logspace(log10(SMI),log10(SMA),SNUM2);
del=floor(bb(2:end)-bb(1:end-1))+1;
[n,x]=histc(P_num',bb);
cc=sqrt(bb(2:end).*bb(1:end-1));%behind
xx=[NUM(:,1)',cc];
yy=[NUM(:,3)'./100,n(1:end-1)'./del/size(P_num',1)];
loglog(xx, yy,'ks-', 'markersize', 5); 
% hold on; loglog(xx,xx.^-1.5,'r-');

indexY0=find(yy==0);%liner fit
IndexXn=[find(xx>NLUp),find(xx<NLDw)];
indexXY=[IndexXn,indexY0];
yy(indexXY)=[];xx(indexXY)=[];
xx=log10(xx);
yy=log10(yy);
K=polyfit(xx,yy,1);slope=K(1);
% xx=min(xx):0.01:max(xx);
% yy=polyval(K,xx);
% hold on;
% loglog(10.^xx,10.^yy,'b','Linewidth',2)
% hold on;
% sumP3=sum( n(1:end-1)'./del/size(P_num',1))+sum(NUM(:,3)'/100);
% axis([0,SMA,10^-5,1.1]);

%  bb=logspace(0,log10(SMA),SNUM1);
%  [n,x]=histc(P_num',bb);
%  loglog(bb(1:end-1), n(1:end-1)'./diff(bb)/size(P_num',1),'oc-', 'markersize', 5);  %在对数坐标轴上 平滑？
%  hold on; loglog(bb(1:end-1), bb(1:end-1).^-1.5,'r--');
%  xx=bb(1:end-1);yy=n(1:end-1)'./diff(bb)/size(P_num',1);
%  indexY0=find(yy==0);%liner fit
%  IndexXn=[find(xx>NLUp),find(xx<NLDw)];
%  indexXY=[IndexXn,indexY0];
%  yy(indexXY)=[];xx(indexXY)=[];
%  xx=log10(xx);
%  yy=log10(yy);
%  K=polyfit(xx,yy,1);slope=K(1);
%  xx=min(xx):0.01:max(xx);
%  yy=polyval(K,xx);
%  hold on;
%  loglog(10.^xx,10.^yy,'r','Linewidth',2)

 
 
 
%  xlabel(['Avalanche Size  ','(\Deltat=',num2str(windw), 'ms','  D=',num2str(D1/100),')'],'FontWeight','demi','FontSize',12);
%  ylabel('Probability','FontWeight','demi','FontSize',12);
%  axis([0,SMA,10^-5,1.1]) 
% % legend('origianl data','reference of -1.5','processed data','fitted');
% legend('reference of -1.5','processed data','fitted line','FontSize',12);

end

