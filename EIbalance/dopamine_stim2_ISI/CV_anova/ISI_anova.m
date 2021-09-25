clear all,close all
%dt=0.1;N=10000;Ne=8000;Ni=N-Ne;f=0.1;
% windw=0.4;
for N_round=1:1:15
    path=['V',num2str(N_round),'\'];
    d1_start=0;d1_end=2;d_d1=0.1;%contrast factor of connection strength
    D1_start=d1_start*100;D1_end=d1_end*100;d_D1=d_d1*100;
    D1_refer=1;

    data0=load([path,'num_parameter_0_',num2str(D1_start),'.log']);
    N=data0(1);f=data0(4);dt=data0(5);life=data0(6);PE=data0(2);
    muEext=data0(7);muIext=data0(8);Tprestim_PT=data0(9);TCamp=data0(10);Tcue_PT=data0(11);
    Ne=round(N*PE/100);Ni=N-Ne;

    size_num=round((D1_end-D1_start)/d_D1)+1;
    VC=zeros(1,size_num);
    recurrent_All=zeros(1,size_num);
    recurrent_E1=zeros(1,size_num);

    n_nostim=round(Tprestim_PT/TCamp);
    n_rest=round((life-Tprestim_PT-Tcue_PT)/TCamp);
    not_count_T=round((Tprestim_PT+Tcue_PT)/TCamp);

    TCamp1=1;
    Tpre=1:Tprestim_PT/TCamp1;Tdur=Tprestim_PT/TCamp1+1:(Tprestim_PT+Tcue_PT)/TCamp1;
    Taft=(Tprestim_PT+Tcue_PT)/TCamp1+1:life/TCamp1;

    num=0;D1=D1_start;
    while D1>=D1_start&&D1<=D1_end
    num=num+1;

    raster=load([path,'rasters_0_',num2str(D1),'.log']);
    % raster(find(raster(:,2)<Tpre(1)|raster(:,2)>Tpre(end)),:)=[];
    raster(find(raster(:,2)<Tdur(1)|raster(:,2)>Tdur(end)),:)=[];
    % raster(find(raster(:,2)<Taft(1)|raster(:,2)>Taft(end)),:)=[];

    Tisi=[];
    for target=1:800
    T_train=raster(find(raster(:,1)==target),2);%ISI
    Tisi=[Tisi;T_train(2:end)-T_train(1:end-1)];
    end

    VC(num)=std(Tisi)/mean(Tisi);%vari coefficient

    data5=load([path,'currents_0_',num2str(D1),'.log']); 
    % data5(find(data5(:,1)<(Tprestim_PT+Tcue_PT)),:)=[];
    % data5(find(data5(:,1)>=(Tprestim_PT+Tcue_PT)),:)=[];
    % recurrent_E1(num)=mean(data5(:,2));
    % recurrent_All(num)=mean(data5(:,2)+data5(:,3)+data5(:,4)+data5(:,5)+data5(:,6)+data5(:,7));

    disp([num2str(D1),'->',num2str(D1_end)]);
    D1=D1+d_D1;
    end

    VC_All(N_round,:)=VC;
    
%     xx=d1_start:d_d1:d1_end;
%     figure();
%     plot(xx,VC,'.-');
%     title('VC','FontWeight','demi','FontSize',12);

end

xx=d1_start:d_d1:d1_end;
figure();
errorbar(xx,mean(VC_All,1),std(VC_All,0,1),'.-');
ylabel('CV','FontWeight','demi','FontSize',12);
xlabel('DA','FontWeight','demi','FontSize',12);
ylim([0.55,0.71]);xlim([0,2]);
%for anvona
[p,table,stats] = anova1(VC_All)
compare=multcompare(stats)
[p,table,stats] = anova1(VC_All)
figure(),boxplot(VC_All)

