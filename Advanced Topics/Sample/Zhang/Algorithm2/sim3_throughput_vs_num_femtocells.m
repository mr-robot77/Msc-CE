clear;
K=200;   % Totall repeat of simulation
nFemto=[4,8,12,16];
nMusr=6;
N=18;
X=1000;
Y=1000;
nFusr=4;

for c=1 : 4
    for j=1 : nFemto(c)+1
        avg_cell_TR(j,c)=0;
    end
end

for c=1 : 4
    for k=1 : K
        macrocell(k) = Macro(nMusr,N,X,Y,nFemto(c),nFusr);
        macrocell(k)=PathGain(macrocell(k));
        macrocell(k)=Update(macrocell(k),1);
        for t=2 : 12
            macrocell(k)=Update(macrocell(k),0);
            
        end
        avg_cell_TR(1,c)=avg_cell_TR(1,c)+ macrocell(k).Cell_TR(t);
        for j=2 : nFemto(c)+1
            avg_cell_TR(j,c)=avg_cell_TR(j,c)+ macrocell(k).Fcel(j-1).Cell_TR(t);
        end
    end
   
    femto_av(c)=0;
    for j=1 : nFemto(c)+1
        avg_cell_TR(j,c)=avg_cell_TR(j,c)/K;
        if j>1
            femto_av(c)=femto_av(c)+ avg_cell_TR(j,c);
        end
    end
    femto_av(c)= femto_av(c)/nFemto(c);
    
    
end
% plot(1:1:c, avg_cell_TR(1,:), 1:1:t,avg_cell_TR(2,:), 1:1:t,avg_cell_TR(3,:),1:1:t,avg_cell_TR(4,:),1:1:t,avg_cell_TR(5,:),1:1:t,avg_cell_TR(6,:),1:1:t,avg_cell_TR(7,:),1:1:t,avg_cell_TR(8,:),1:1:t,avg_cell_TR(9,:));
% xlabel('Time(Iteration)')
% ylabel('Throughput(bps/HZ)')
% legend('Macrocell','Femtocell 1','Femtocell 2','Femtocell 3','Femtocell 4','Femtocell 5','Femtocell 6','Femtocell 7','Femtocell 8')
% grid
% figure

plot(1:1:c, avg_cell_TR(1,:), 1:1:c,femto_av(:));
xlabel('Time(Iteration)')
ylabel('Throughput(bps/HZ)')
legend('Macrocell','Femtocells Average ');
grid
figure
