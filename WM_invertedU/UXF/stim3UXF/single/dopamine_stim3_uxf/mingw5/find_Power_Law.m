function [expect] = find_Power_Law(raster,windw,life,Ne,pr,f,Jc)
p=5;
%%%%%%%%%%%%%%%%%%%for LFP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% th_n=2;N_lfp=100;

% raster(find(raster(:,1)>Ne),:)=[];%ֻȡ�˷���Eϸ��
raster(find(raster(:,1)>=Ne*p*f),:)=[];%ֻȡ�˷���Eϸ����ѡ����ϸ��
% raster(find(raster(:,1)>Ne*f),:)=[];%ֻȡ�˷���E1ϸ����ѡ����ϸ��

index_r=find(mod(raster(:,1),round(1/pr))~=0);%��Ŀ�꼯��ȡpr*100*%���ֹ۲�
raster(index_r,:)=[];

%%%%%%%%%%%%%%%%%%%for LFP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raster(:,1)=floor(raster(:,1)./N_lfp);raster(:,3:4)=[];
% raster(:,2)=floor(raster(:,2)./windw); 
% raster(find(raster(:,2)==max(raster(:,2))),:)=[];%ȥβ
% ss=Ne*p*f/N_lfp;ras_LFP=zeros(ss,max(raster(:,2))+1);
% for i=1:numel(raster(:,2))
%     ras_LFP(raster(i,1)+1,raster(i,2)+1)=ras_LFP(raster(i,1)+1,raster(i,2)+1)+1;
% end
% r_lfp=(ras_LFP>=th_n);
% LFP=sum(r_lfp,1);
% fire_num=[0,LFP,0];


raster_t=raster(:,2);
% raster_t=raster(1:round(end/10),2);
life=floor(life/windw)*windw;
central_window=(windw/2):windw:life;
fire_num=hist(raster_t,central_window);%��ÿ��������Ƶ��
fire_num=[0,fire_num,0];%ȷ����βΪ0


index_0=find(fire_num==0);   
index_inter=index_0(:,2:size(index_0,2))-index_0(:,1:size(index_0,2)-1);%��0Ԫ�ؼ���
index_inter(size(index_inter,2))=index_inter(size(index_inter,2))+1;%���һ���ذ����0����
P_cell=mat2cell(fire_num,1,index_inter);

P_nc=cell(size(P_cell));
for i=1:size(index_inter,2) %���
   P_nc(i)={sum(cat(2,P_cell{1,i},0))};   
end
   P_num=cell2mat(P_nc);
   P_num(P_num==0)=[];%ȥ0

 m=0;n_max=max(P_num);  
for k=1:size(index_0,2)-1
    if fire_num(index_0(k)+1)~=0
         m=m+1;
%           expect(m)=(fire_num(index_0(k)+2)/fire_num(index_0(k)+1))*(n_max-1)/(n_max-fire_num(index_0(k)+1));
          expect(m)=round(fire_num(index_0(k)+2)/fire_num(index_0(k)+1))*(n_max-1)/(n_max-fire_num(index_0(k)+1));
%          expect(m)=(fire_num(index_0(k)+2)/fire_num(index_0(k)+1));
         
    end
end
expect=mean(expect);

 bb=logspace(0,log10(max(P_num')),30);
 [n,x]=histc(P_num',bb);
 loglog(bb(1:end-1), n(1:end-1)'./diff(bb)/size(P_num',1),'o-', 'markersize', 5);  %�ڶ����������� ƽ����
 hold on; loglog(bb(1:end-1), bb(1:end-1).^-1.5,'r--');
 hold on; loglog(bb(1:end-1), bb(1:end-1).^-2.1,'b--');hold off;
 xlabel(['Claster Size  ','(Window=',num2str(windw), 'ms','  Jc=',num2str(Jc/100),')'],'FontWeight','demi','FontSize',12);
 ylabel('Probability','FontWeight','demi','FontSize',12);
   
end

%%%test%%%%% 
% life=23;
% window=4;
% c=randi(20,20,1);
% d=c'
% life=floor(life/window)*window
% core_point=(window/2):window:life
% n=hist(c,core_point)

% clear all;
% c=randi(5,1,20);xx=ceil(rand(1,20)*20);c(xx)=0;
% fire_num=c
% fire_num=[0,fire_num,0];%ȷ����βΪ0
% index_0=find(fire_num==0);
% index_inter=index_0(:,2:size(index_0,2))-index_0(:,1:size(index_0,2)-1);
% index_inter(size(index_inter,2))=index_inter(size(index_inter,2))+1;
% P_cell=mat2cell(fire_num,1,index_inter);
% P_nc=cell(size(P_cell));
% for i=1:size(index_inter,2) 
%    P_nc(i)={sum(cat(2,P_cell{1,i},0))};   
% end
%    P_num=cell2mat(P_nc)
%    P_num(P_num==0)=[]
% 
%  figure(); x=1:10:1000;y=0.1:10:10;loglog(x,y,'*');
