clear all
%dt=0.1;N=10000;Ne=8000;Ni=N-Ne;f=0.1;
% windw=0.4;
d1_start=0;d1_end=2;d_d1=0.1;%dopamine concentration
D1_start=d1_start*100;D1_end=d1_end*100;d_D1=d_d1*100;

size_num=round((D1_end-D1_start)/d_D1)+1;

figure();
color=['b','k','g'];
for nn=1:2%1:1,1:2;

minStim=zeros(1,size_num);
num=0;D1=D1_start;
while D1>=D1_start&&D1<=D1_end
num=num+1;
  
data0=load(['.\dopamine_stim2_min_stim_strength',num2str(nn),'\mingw5\num_parameter_0_',num2str(D1),'.log']);
N=data0(1);f=data0(4);dt=data0(5);life=data0(6);PE=data0(2);
muEext=data0(7);muIext=data0(8);Tprestim_PT=data0(9);TCamp=data0(10);Tcue_PT=data0(11);
Ne=round(N*PE/100);Ni=N-Ne;

minStim(num)=Tcue_PT;%min stim time
D1=D1+d_D1;
end
% minStim(15)=45;
sens=1./minStim*100;%sensitivity

xx=d1_start:d_d1:d1_end;
%%%%%%%%%%%%%%%% min stim time,sensitivity%%%%%%%%%

subplot(1,2,1);
plot(xx,minStim,[color(nn),'.-'],'LineWidth',1.5,'MarkerSize',16);
% title('min stimulation time','FontWeight','demi','FontSize',12);
xlabel('Dopamine','FontWeight','demi','FontSize',12); 
ylabel('T_{min}','FontWeight','demi','FontSize',12);hold on;

subplot(1,2,2);
plot(xx,sens,[color(nn),'.-'],'LineWidth',1.5,'MarkerSize',16);
% title('sensitivity for stimulation','FontWeight','demi','FontSize',12);
xlabel('Dopamine','FontWeight','demi','FontSize',12); 
ylabel('[%]','FontWeight','demi','FontSize',12); hold on;
 ylim([0,30]);
end
subplot(1,2,1);legend('strength1','strength2');
subplot(1,2,2);legend('strength1','strength2');
