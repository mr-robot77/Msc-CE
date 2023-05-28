clear;
K=50;   % Totall repeat of simulation
nFemto=8;
nMusr=6;
N=18;
X=1000;
Y=1000;
nFusr=4;
macro_Pmx=[0.01,0.01,0.01,0.01];
% macro_Tsnr=[9.55,9.55,9.55,9.55];
macro_minRate=[512000,1024000,1536000,2048000];
macro_QAM=[1,2,3,4];
for k=1 : K
        macrocell(k) = Macro(nMusr,N,X,Y,nFemto,nFusr);
end
T=7;
for c=1 : 4
    for j=1 : nFemto+1
        cell_TR(j,c)=0;
    end
end
for c=1 : 4
%     macrocell(k)=Reset_It(macrocell(k));
    for k=1 : K
        macrocell(k)=Reset_It(macrocell(k));
%         macrocell(k) = Macro(nMusr,N,X,Y,nFemto,nFusr);
        macrocell(k)=PathGain(macrocell(k));
        macrocell(k)=Set_conf2(macrocell(k),macro_minRate(c),macro_QAM(2));
        macrocell(k)=Update(macrocell(k),1);
        
        for t=2 : T
            macrocell(k)=Update(macrocell(k),0);
        end
        
        cell_TR(1,c)=cell_TR(1,c)+ macrocell(k).Cell_TR(T);
        for j=2 : nFemto+1
            cell_TR(j,c)=cell_TR(j,c)+ macrocell(k).Fcel(j-1).Cell_TR(T);
        end
    end
    femto_av(c)=0;
    for j=1 : nFemto+1
        cell_TR(j,c)=cell_TR(j,c)/K;
        if j>1
            femto_av(c)=femto_av(c)+ cell_TR(j,c);
        end
    end
end

plot(1:1:c, cell_TR(1,:), 1:1:c,femto_av(:));
xlabel('Time(Iteration)')
ylabel('Throughput(bps/HZ)')
legend('Macrocell','Femtocells Average ');
grid


%     plot(1:1:c, cell_TR(1,:), 1:1:c,cell_TR(2,:), 1:1:c,cell_TR(3,:),1:1:c,cell_TR(4,:),1:1:c,cell_TR(5,:),1:1:c,cell_TR(6,:),1:1:c,cell_TR(7,:),1:1:c,cell_TR(8,:),1:1:c,cell_TR(9,:));
%     xlabel('Macrocell Users Target SINR')
%     ylabel('Throughput')
%     legend('Macrocell','Femtocell 1','Femtocell 2','Femtocell 3','Femtocell 4','Femtocell 5','Femtocell 6','Femtocell 7','Femtocell 8');
%     grid
