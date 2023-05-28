clear;
K=1;   % Totall repeat of simulation
nFemto=8;
nMusr=4;
N=8;
X=1000;
Y=1000;
nFusr=4;
macro_Pmx=[0.01,0.01,0.01,0.01];
% macro_Tsnr=[9.55,9.55,9.55,9.55];
macro_Tsnr=[9.55,45.11,179.85,694.17];
for k=1 : K
        macrocell(k) = Macro(nMusr,N,X,Y,nFemto,nFusr);
end
for c=1 : 4
    for j=1 : nFemto+1
        cell_TR(j,c)=0;
    end
end
for c=1 : 4
    macrocell(k)=Reset_It(macrocell(k));
    for k=1 : K
        macrocell(k)=PathGain(macrocell(k));
        macrocell(k)=Set_conf(macrocell(k),macro_Pmx(c),macro_Tsnr(c));
        macrocell(k)=Update(macrocell(k),1);
        
        for t=2 : 20
            macrocell(k)=Update(macrocell(k),0);
        end
        
        cell_TR(1,c)=cell_TR(1,c)+ macrocell(k).Cell_TR(20);
        for j=2 : nFemto+1
            cell_TR(j,c)=cell_TR(j,c)+ macrocell(k).Fcel(j-1).Cell_TR(20);
        end
    end
    for j=1 : nFemto+1
        cell_TR(j,c)=cell_TR(j,c)/K;
    end
end

    plot(1:1:c, cell_TR(1,:), 1:1:c,cell_TR(2,:), 1:1:c,cell_TR(3,:),1:1:c,cell_TR(4,:),1:1:c,cell_TR(5,:),1:1:c,cell_TR(6,:),1:1:c,cell_TR(7,:),1:1:c,cell_TR(8,:),1:1:c,cell_TR(9,:));
    xlabel('Macrocell Users Target SINR')
    ylabel('Throughput')
    legend('Macrocell','Femtocell 1','Femtocell 2','Femtocell 3','Femtocell 4','Femtocell 5','Femtocell 6','Femtocell 7','Femtocell 8');
    grid
