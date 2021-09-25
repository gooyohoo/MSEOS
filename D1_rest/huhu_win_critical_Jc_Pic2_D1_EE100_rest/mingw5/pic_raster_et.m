function [] = pic_raster_et( data1,data2,data4,data3,data5,dt,Ne,f,D1 )
%PIC_RASTER_ET Summary of this function goes here
%   Detailed explanation goes here

subplot(2,2,1);
plot(data1(:,1),data1(:,2),'r',data1(:,1),data1(:,7),'b');
title('fire rate of the Exc1&Inh  ','FontWeight','demi','FontSize',12);
xlabel(['TIME  ( \Delta','t=',num2str(dt), 'ms)'],'FontWeight','demi','FontSize',12); 
ylabel('fire rate[Hz]','FontWeight','demi','FontSize',12);
legend(['Fir_E ','Fir_I:']);

subplot(2,2,2);
plot(data2(:,1),data2(:,2),'r');hold on;
plot(data3(:,1),data3(:,2),'b');hold on;
title('facilitation: <u(t)> & depression: <x(t)>','FontWeight','demi','FontSize',12);
xlabel(['TIME  ( \Delta','t=',num2str(dt), 'ms)'],'FontWeight','demi','FontSize',12);
ylabel('Value','FontWeight','demi','FontSize',12); 
legend('u: Fac.','x: Dep.');


subplot(2,2,3);
plot(data5(:,1),data5(:,2),'r',data5(:,1),data5(:,7),'b');
title('Current of a single neuron','FontWeight','demi','FontSize',12);
xlabel(['TIME  ( \Delta','t=',num2str(dt), 'ms)'],'FontWeight','demi','FontSize',12);
ylabel('Value','FontWeight','demi','FontSize',12); 
legend('I_e','I_i');   

subplot(2,2,4);
data4(:,1)=data4(:,1)/10;
siz=size(data4,1);
data4E=data4.*[ones(siz,3),ones(siz,1)];
index=find(data4E(:,3)~=1);data4E(index,:)=[];
data4E1=data4E;data4E2=data4E;
index_E1=find(data4E1(:,4)==0);data4E1(index_E1,:)=[];
index_E21=find(data4E2(:,4)==1);data4E2(index_E21,:)=[];
% index_E22=find(data4E2(:,1)>=Ne*f*0.1*2);
index_E22=find(data4E2(:,1)>=Ne*f*0.1*5);
data4E2(index_E22,:)=[];data4E2(:,1)=data4E2(:,1);%-f*Ne*0.1;
plot(data4E1(:,2)',data4E1(:,1)','r.',data4E2(:,2)',data4E2(:,1)','b.');
title('Raster activity of E1&E2 ','FontWeight','demi','FontSize',12);
xlabel(['D1=',num2str(D1)],'FontWeight','demi','FontSize',12); 
ylabel('The neuron index','FontWeight','demi','FontSize',12);
ylim([0,400]);
legend('E1','E_r');
end

