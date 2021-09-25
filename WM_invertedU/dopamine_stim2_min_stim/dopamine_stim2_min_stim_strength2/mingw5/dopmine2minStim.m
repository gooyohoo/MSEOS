clear all
%dt=0.1;N=10000;Ne=8000;Ni=N-Ne;f=0.1;
% windw=0.4;
d1_start=0;d1_end=2;d_d1=0.1;%dopamine concentration
D1_start=d1_start*100;D1_end=d1_end*100;d_D1=d_d1*100;

size_num=round((D1_end-D1_start)/d_D1)+1;
minStim=zeros(1,size_num);

num=0;D1=D1_start;
while D1>=D1_start&&D1<=D1_end
num=num+1;
  
data0=load(['num_parameter_0_',num2str(D1),'.log']);
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
figure();
subplot(1,2,1);
plot(xx,minStim,'s-');
title('min stim time','FontWeight','demi','FontSize',12);
xlabel('dopamine D1','FontWeight','demi','FontSize',12); 
ylabel('t[ms]','FontWeight','demi','FontSize',12);
grid;
subplot(1,2,2);
plot(xx,sens,'o-');
title('sensitivity for stim','FontWeight','demi','FontSize',12);
xlabel('dopamine D1','FontWeight','demi','FontSize',12); 
ylabel('%','FontWeight','demi','FontSize',12);
% ylim([0.6,1.1]);
grid;
