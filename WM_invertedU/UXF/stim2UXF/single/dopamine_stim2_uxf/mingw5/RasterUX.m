clear all
D1=100;%contrast factor of connection strength
windw=0.3;% time bin

data0=load(['num_parameter_0_',num2str(D1),'.log']);
N=data0(1);PE=data0(2);D1=data0(3);f=data0(4);dt=data0(5);life=data0(6);
muEext=data0(7);muIext=data0(8);Tprestim_PT=data0(9);TCamp=data0(10);
Ne=round(N*PE/100);Ni=N-Ne;

data2=load(['stp_u_0_',num2str(D1),'.log']);
data3=load(['stp_x_0_',num2str(D1),'.log']);
raster=load(['rasters_0_',num2str(D1),'.log']);
index_data4=find(mod(raster(:,1),10)~=0);
data4=raster;data4(index_data4,:)=[];

data4(:,1)=data4(:,1)/10;
siz=size(data4,1);
data4E=data4.*[ones(siz,3),ones(siz,1)];
index=find(data4E(:,3)~=1);data4E(index,:)=[];
data4E1=data4E;data4E2=data4E;
index_E1=find(data4E1(:,4)==0);data4E1(index_E1,:)=[];%接受刺激的
index_E21=find(data4E2(:,4)==1);data4E2(index_E21,:)=[];%不接受刺激的
index_E22=find(data4E2(:,1)>=Ne*f*0.1*2);data4E2(index_E22,:)=[];


figure();
rectangle('Position',[Tprestim_PT,0,200,80],'EdgeColor','w','FaceColor',[0.8,0.8,0.8])
hold on;

t1=data2(:,1);%ux
u1=data2(:,2);u2=data2(:,3)+1;
x1=data3(:,2);x2=data3(:,3)+1;
t21=data4E1(:,2)';t22=data4E2(:,2)';
sp1=data4E1(:,1)';sp2=data4E2(:,1)';
[AX,H1,H2]=plotyy([t21,t22],[sp1,sp2],t1,[u2,u1,x1,x2]*0.8,'plot');

set(H1,'LineStyle','o','Color','k','Markersize',3,'markerfacecolor','k');
set(H2(1),'Color','b','Linewidth',2);
set(H2(2),'Color','b','Linewidth',2);
set(H2(3),'Color','r','Linewidth',2);
set(H2(4),'Color','r','Linewidth',2);

set(AX(1),'XColor','k','YColor','k');
set(AX(2),'XColor','k','YColor','k')

set(AX(2),'YTick',(0:0.5:2)*0.8,'YTickLabel',{'0','0.5','1（0）','0.5','1'})
set(AX(1),'YTick',0:80:160,'YTickLabel',{'0','80','160'})
set(AX,'XTick',0:1000:4000,'XTickLabel',{'','1000','2000','3000',''});

HH1=get(AX(1),'Ylabel');
set(HH1,'String','Cell Index','FontWeight','demi','FontSize',12);
HH2=get(AX(1),'Xlabel');
set(HH2,'String','t[ms]','FontWeight','demi','FontSize',12);

legend([H2(1),H2(3)],{'x','u'},'Color','w');

