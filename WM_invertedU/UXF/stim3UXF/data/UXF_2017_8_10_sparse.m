clear all
%dt=0.1;N=10000;Ne=8000;Ni=N-Ne;f=0.1;
% windw=0.4;

Trial_end=9;Trial_start=1;

d1_start=0;d1_end=2;d_d1=0.1;%contrast factor of connection strength
D1_start=d1_start*100;D1_end=d1_end*100;d_D1=d_d1*100;
D1_refer=1;
size_num=round((D1_end-D1_start)/d_D1)+1;

data0=load(['./dopamine_stim3_uxf1/mingw5/','num_parameter_0_',num2str(D1_start),'.log']);
N=data0(1);f=data0(4);dt=data0(5);life=data0(6);PE=data0(2);
muEext=data0(7);muIext=data0(8);Tprestim_PT=data0(9);TCamp=data0(10);Tcue_PT=data0(11);
Ne=round(N*PE/100);Ni=N-Ne;

M2=(d1_end-d1_start)/d_d1;
M1=Trial_end-Trial_start;
u_stim_end=zeros(M1,M2);
x_stim_end=zeros(M1,M2);
u_m=zeros(M1,M2);
x_m=zeros(M1,M2);

Fe1_stim_end=zeros(M1,M2);
Fe1_stim_dur=zeros(M1,M2);
Fe1_m=zeros(M1,M2);
Fe1_m2=zeros(M1,M2);
n_nostim=round(Tprestim_PT/TCamp);
n_rest=round((life-Tprestim_PT-Tcue_PT)/TCamp);
not_count_T=round((Tprestim_PT+Tcue_PT)/TCamp);

for num1=Trial_start:Trial_end
     
    path=['./dopamine_stim3_uxf',num2str(num1),'/mingw5/'];
    
    num2=0;D1=D1_start;
    while D1>=D1_start&&D1<=D1_end
        num2=num2+1;

        data2=load([path,'stp_u_0_',num2str(D1),'.log']);
        u_stim_end(num1,num2)=mean(mean(data2(not_count_T,1:800)));

        data3=load([path,'stp_x_0_',num2str(D1),'.log']);
        x_stim_end(num1,num2)=mean(mean(data3(not_count_T,1:800)));

        data4=load([path,'rates_pops_0_',num2str(D1),'.log']);
        Fe1_stim_end(num1,num2)=mean(data4(not_count_T,2));
        Fe1_stim_dur(num1,num2)=mean(data4(n_nostim:not_count_T,2));
        
        Fe1_m(num1,num2)=mean(data4(not_count_T:end,2));
        
        Fe1_m2(num1,num2)=mean(data4(not_count_T+120:end,2));
        
        u_m(num1,num2)=mean(mean(data2(not_count_T+100:end,1:800)));
        x_m(num1,num2)=mean(mean(data3(not_count_T+100:end,1:800)));


        D1=D1+d_D1;
    end
    disp(num1);

end

ux_stim_end=u_stim_end.*x_stim_end;

%Fe1=Fe1_stim_dur;
%Fe1=Fe1_stim_end;
Fe1=Fe1_m;
%Fe1=Fe1_m2;

uxf=ux_stim_end.*Fe1;

        
% ux_m=u_m.*x_m;
% uxf=ux_m.*Fe1;

xx=d1_start:d_d1:d1_end;
OptiDA=1/d_d1+1;

MFe1=mean(Fe1,1);SFe1=(std(Fe1,0,1));
Mu_stim_end=mean(u_stim_end,1);Su_stim_end=std(u_stim_end,0,1);
Mx_stim_end=mean(x_stim_end,1);Sx_stim_end=std(x_stim_end,0,1);
Mux_stim_end=mean(ux_stim_end,1);Sux_stim_end=std(ux_stim_end,0,1);
Muxf=mean(uxf,1);Suxf=std(uxf,0,1);

%% 
% temp1=[1,2,3,4,5,6,7,8,9,11,10,12,13,14,15,16,17,18,19,20,21];
% temp2=[1,2,3,4,5,6,7,8,9,10,12,11,13,14,15,16,17,18,19,20,21];
% Mu_stim_end=Mu_stim_end(temp1);Su_stim_end=Su_stim_end(temp1);
% Mx_stim_end=Mx_stim_end(temp1);Sx_stim_end=Sx_stim_end(temp1);
% Mux_stim_end=Mux_stim_end(temp1);Sux_stim_end=Sux_stim_end(temp1);
% Muxf=Muxf(temp2);Suxf=Suxf(temp2);
% MFe1=MFe1(temp2);SFe1=SFe1(temp2);
%%

figure();
subplot(1,5,1);
errorbar(xx,MFe1,SFe1,'-','color','k');
xlabel('DA','FontWeight','demi','FontSize',12); 
ylabel('f[Hz]','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0,MFe1(OptiDA)],'line','--','color','k');hold off;
ylim([0,16]);xlim([-0.2,2.2]);
subplot(1,5,2);
errorbar(xx,Mu_stim_end,Su_stim_end,'-','color','k');
xlabel('DA','FontWeight','demi','FontSize',12); 
ylabel('u','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0,Mu_stim_end(OptiDA)],'line','--','color','k');hold off;
ylim([0,1]);xlim([-0.2,2.2]);
subplot(1,5,3);
errorbar(xx,Mx_stim_end,Sx_stim_end,'-','color','k');
xlabel('DA','FontWeight','demi','FontSize',12); 
ylabel('x','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0,Mx_stim_end(OptiDA)],'line','--','color','k');hold off;
ylim([0,0.8]);xlim([-0.2,2.2]);
subplot(1,5,4);
errorbar(xx,Mux_stim_end,Sux_stim_end,'-','color','k');
xlabel('DA','FontWeight','demi','FontSize',12); 
ylabel('ux','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0,Mux_stim_end(OptiDA)],'line','--','color','k');hold off;
ylim([0,0.4]);xlim([-0.2,2.2]);
subplot(1,5,5);
errorbar(xx,Muxf,Suxf,'-','color','k');
xlabel('DA','FontWeight','demi','FontSize',12); 
ylabel('uxf','FontWeight','demi','FontSize',12);
hold on;line([D1_refer,D1_refer],[0,Muxf(OptiDA)],'line','--','color','k');hold off;
ylim([0,4]);xlim([-0.2,2.2]);

% figure();
% subplot(1,5,1);
% plot(xx,mean(Fe1_m,1),'-');
% xlabel('DA','FontWeight','demi','FontSize',12); 
% ylabel('f[Hz]','FontWeight','demi','FontSize',12);
% hold on;line([D1_refer,D1_refer],[0,mean(Fe1_m(:,OptiDA))],'line','--');hold off;
% ylim([0,20]);
% subplot(1,5,2);
% plot(xx,mean(u_stim_end,1),'-');
% xlabel('DA','FontWeight','demi','FontSize',12); 
% ylabel('u','FontWeight','demi','FontSize',12);
% hold on;line([D1_refer,D1_refer],[0,mean(u_stim_end(:,OptiDA))],'line','--');hold off;
% ylim([0,1]);
% subplot(1,5,3);
% plot(xx,mean(x_stim_end,1),'-');
% xlabel('DA','FontWeight','demi','FontSize',12); 
% ylabel('x','FontWeight','demi','FontSize',12);
% hold on;line([D1_refer,D1_refer],[0,mean(x_stim_end(:,OptiDA))],'line','--');hold off;
% ylim([0,0.8]);
% subplot(1,5,4);
% plot(xx,mean(ux_stim_end,1),'-');
% xlabel('DA','FontWeight','demi','FontSize',12); 
% ylabel('ux','FontWeight','demi','FontSize',12);
% hold on;line([D1_refer,D1_refer],[0.1,mean(ux_stim_end(:,OptiDA))],'line','--');hold off;
% ylim([0.1,0.4]);
% subplot(1,5,5);
% plot(xx,mean(uxf,1),'-');
% xlabel('DA','FontWeight','demi','FontSize',12); 
% ylabel('uxf','FontWeight','demi','FontSize',12);
% hold on;line([D1_refer,D1_refer],[0.1,mean(uxf(:,OptiDA))],'line','--');hold off;
% ylim([0,4]);


