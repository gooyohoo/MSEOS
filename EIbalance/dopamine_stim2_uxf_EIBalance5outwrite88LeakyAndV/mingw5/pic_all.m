function [] = pic_all( data1,data2,data4,data3,data5,dt,Ne,f,D1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p=5;
figure();
for pi=1:1:p
subplot(2,3,pi);
plot(data1(:,1),data1(:,pi+1),'r');
title(['fire rate of the Exc',num2str(pi)],'FontWeight','demi','FontSize',12);
xlabel(['D1=',num2str(D1)],'FontWeight','demi','FontSize',12); 
ylabel('fire rate[Hz]','FontWeight','demi','FontSize',12);
legend('Fir_E',num2str(pi));
end
subplot(2,3,6);
plot(data1(:,1),data1(:,7),'b');
title('fire rate of the Inh  ','FontWeight','demi','FontSize',12);
xlabel(['D1=',num2str(D1)],'FontWeight','demi','FontSize',12); 
ylabel('fire rate[Hz]','FontWeight','demi','FontSize',12);
legend('Fir_I');

figure();
for pi=1:1:p
subplot(2,3,pi);
plot(data2(:,1),data2(:,pi+1),'r');hold on;
plot(data3(:,1),data3(:,pi+1),'b');hold on;
title([' <u(t)>&<x(t)>',num2str(pi)],'FontWeight','demi','FontSize',12);
xlabel(['D1=',num2str(D1)],'FontWeight','demi','FontSize',12); 
ylabel('Value','FontWeight','demi','FontSize',12); 
legend('u: Fac.','x: Dep.');
end

figure();
for pi=1:1:p
subplot(2,3,pi);
plot(data5(:,1),data5(:,pi+1),'r');
title('recurrent of pop per neuron','FontWeight','demi','FontSize',12);
xlabel(['D1=',num2str(D1)],'FontWeight','demi','FontSize',12);
ylabel('Value','FontWeight','demi','FontSize',12); 
legend(['I_e',num2str(pi)]);   
end
subplot(2,3,6);
plot(data5(:,1),data5(:,7),'b');
title('recurrent of pop per neuron','FontWeight','demi','FontSize',12);
xlabel(['D1=',num2str(D1)],'FontWeight','demi','FontSize',12);
ylabel('Value','FontWeight','demi','FontSize',12); 
legend('I_i');  


figure();
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
legend('E1','E2');
end




